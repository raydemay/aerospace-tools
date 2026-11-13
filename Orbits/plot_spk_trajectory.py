"""
SPICE SPK 3D Trajectory Visualizer
---------------------------------
Reads a SPICE SPK ephemeris kernel (.bsp / .spk) using SpiceyPy, computes
the trajectory of a spacecraft or celestial body relative to the Sun, and
renders a 3D visualization using both Matplotlib and Plotly.

Features:
- Auto-detects NAIF body IDs and coverage intervals within the SPK file.
- Heliocentric J2000 Ecliptic frame ('ECLIPJ2000') positioning.
- Plots the spacecraft/body trajectory, Earth's orbital path, and the Sun.
- Outputs an interactive Plotly HTML file (zoom, pan, rotate) and a Matplotlib 3D view.
"""

import argparse
import math
import os
import sys
from typing import Any

# Astronomical Unit in kilometers (IAU standard)
KM_PER_AU = 149597870.7


def check_dependencies():
    """Verify that required 3rd-party libraries are installed."""
    missing = []
    for pkg in ["spiceypy", "numpy", "matplotlib", "plotly"]:
        try:
            __import__(pkg)
        except ImportError:
            missing.append(pkg)
    if missing:
        print(f"\n[ERROR] Missing required packages: {', '.join(missing)}")
        print("Please install them using pip:")
        print("    pip install spiceypy matplotlib plotly numpy\n")
        sys.exit(1)


def load_spice_kernels(spk_path, extra_kernels=None):
    """
    Load the target SPK file and any additional companion kernels (LSK, PCK, BSP).
    """
    import spiceypy as spice

    spice.kclear()

    # Load extra helper kernels first (e.g. leapseconds .tls or planetary ephemeris .bsp)
    if extra_kernels:
        for k in extra_kernels:
            if os.path.isfile(k):
                spice.furnsh(k)
                print(f"[SPICE] Loaded kernel: {k}")

    # Auto-load any .tls (leapseconds) kernel in the same directory if present
    spk_dir = os.path.dirname(os.path.abspath(spk_path))
    for f in os.listdir(spk_dir):
        if f.lower().endswith(".tls"):
            tls_path = os.path.join(spk_dir, f)
            spice.furnsh(tls_path)
            print(f"[SPICE] Auto-loaded leapseconds kernel: {tls_path}")
            break

    # Load main SPK file
    spice.furnsh(spk_path)
    print(f"[SPICE] Loaded target SPK kernel: {spk_path}")


def inspect_spk_bodies(spk_path):
    """
    Queries the SPK file for all contained NAIF IDs and their coverage windows.
    Returns a list of dicts with body info: [{'id': id, 'name': name, 'start_et': et0, 'end_et': et1}]
    """
    import spiceypy as spice

    # spkobj returns a SPICE cell of integer IDs
    body_ids = spice.spkobj(spk_path)
    bodies = []

    for bid in body_ids:
        # spkcov returns coverage intervals in ET seconds past J2000
        cover = spice.spkcov(spk_path, bid)
        if len(cover) >= 2:
            start_et = cover[0]
            end_et = cover[1]
            try:
                name = spice.bodc2n(bid)
            except Exception:
                name = f"NAIF_ID_{bid}"

            bodies.append({
                "id": bid,
                "name": name,
                "start_et": start_et,
                "end_et": end_et,
            })

    return bodies


def format_et(et):
    """Formats ephemeris time (ET) as human-readable date string."""
    import spiceypy as spice

    try:
        # If an LSK kernel is furnished, timout/et2utc provides exact calendar string
        return spice.timout(et, "YYYY-MON-DD HR:MN:SC (TDB)")
    except Exception:
        # Approximate calendar year if no LSK kernel is loaded
        # 1 Julian year ~ 31557600 TDB seconds past J2000 (2000-01-01 12:00:00)
        year = 2000.0 + et / 31557600.0
        return f"ET: {et:.1f} s (~{year:.2f})"


def sample_trajectory(target_id, start_et, end_et, num_points=500, observer="SUN", frame="ECLIPJ2000"):
    """
    Computes [x, y, z] trajectory positions and dates for target body.
    Returns: (positions array of shape (N, 3) in km, array of et times)
    """
    import numpy as np
    import spiceypy as spice

    ets = np.linspace(start_et, end_et, num_points)
    positions = []
    valid_ets = []

    for et in ets:
        try:
            # spkpos returns (pos, ltime) in km
            pos, _ = spice.spkpos(str(target_id), et, frame, "NONE", observer)
            positions.append(pos)
            valid_ets.append(et)
        except Exception:
            continue

    if not positions:
        raise RuntimeError(
            f"Could not evaluate position for target '{target_id}' relative to '{observer}' "
            f"in frame '{frame}'. Ensure the SPK covers this target."
        )

    return np.array(positions), np.array(valid_ets)


def sample_earth_orbit(start_et, end_et, num_points=200, observer="SUN", frame="ECLIPJ2000"):
    """
    Evaluates Earth's orbital path around the Sun.
    Tries SPICE spkpos first; falls back to an analytical Keplerian 1 AU orbit if
    planetary ephemeris (e.g. de440s.bsp) is not furnished.
    """
    import numpy as np
    import spiceypy as spice

    ets = np.linspace(start_et, end_et, num_points)
    positions = []

    # Attempt to query SPICE for Earth (NAIF ID 399) or Earth-Moon Barycenter (3)
    spice_success = False
    for target in ["399", "3", "EARTH"]:
        try:
            test_pos, _ = spice.spkpos(target, start_et, frame, "NONE", observer)
            # If successful, sample across the interval
            for et in ets:
                pos, _ = spice.spkpos(target, et, frame, "NONE", observer)
                positions.append(pos)
            spice_success = True
            break
        except Exception:
            continue

    if spice_success:
        return np.array(positions)

    # Fallback: Analytical Earth orbit in ECLIPJ2000 (1 AU circle on XY plane)
    print("[INFO] Planetary SPK not loaded; using analytical 1 AU Earth orbit for reference.")
    theta = np.linspace(0, 2 * math.pi, num_points)
    r_earth = KM_PER_AU  # ~149.6 million km
    x = r_earth * np.cos(theta)
    y = r_earth * np.sin(theta)
    z = np.zeros_like(theta)
    return np.column_stack((x, y, z))


def plot_matplotlib(body_name, traj_km, earth_km, ets, unit="AU", save_png=None):
    """
    Creates a 3D trajectory plot using Matplotlib.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    scale = KM_PER_AU if unit.upper() == "AU" else 1.0
    u_label = "AU" if unit.upper() == "AU" else "km"

    traj = traj_km / scale
    earth = earth_km / scale

    fig = plt.figure(figsize=(10, 8), dpi=100)
    ax: Any = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("#0b0e14")
    fig.patch.set_facecolor("#0b0e14")

    # 1. Plot Sun at the origin
    ax.scatter([0], [0], [0], color="#ffcc00", s=180, edgecolors="orange", label="Sun", zorder=10)

    # 2. Plot Earth's orbital path
    ax.plot(earth[:, 0], earth[:, 1], earth[:, 2], color="#29b6f6", linestyle="--", linewidth=1.2, alpha=0.7, label="Earth Orbit")

    # 3. Plot Target Trajectory
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], color="#ff5252", linewidth=2.0, label=f"{body_name} Trajectory")

    # 4. Mark Start and End points
    ax.scatter([traj[0, 0]], [traj[0, 1]], [traj[0, 2]], color="#69f0ae", s=60, marker="o", label="Start", zorder=9)
    ax.scatter([traj[-1, 0]], [traj[-1, 1]], [traj[-1, 2]], color="#e040fb", s=60, marker="^", label="End", zorder=9)

    # Axis styling
    ax.set_xlabel(f"X ({u_label})", color="white")
    ax.set_ylabel(f"Y ({u_label})", color="white")
    ax.set_zlabel(f"Z ({u_label})", color="white")
    ax.tick_params(colors="white")

    start_date = format_et(ets[0])
    end_date = format_et(ets[-1])
    ax.set_title(
        f"{body_name} Heliocentric Trajectory (ECLIPJ2000)\n{start_date} to {end_date}",
        color="white",
        fontsize=12,
        pad=15,
    )

    # Transparent pane backgrounds
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor("#37474f")

    ax.grid(True, color="#37474f", linestyle=":", alpha=0.6)
    legend = ax.legend(facecolor="#1c2430", edgecolor="#37474f", labelcolor="white", loc="upper right")
    legend.get_frame().set_alpha(0.85)

    # Set equal aspect ratio
    max_range = np.array([
        traj[:, 0].max() - traj[:, 0].min(),
        traj[:, 1].max() - traj[:, 1].min(),
        traj[:, 2].max() - traj[:, 2].min(),
        2.5,  # ensure at least 2.5 AU range so Earth orbit is visible
    ]).max() / 2.0

    mid_x = traj[:, 0].mean()
    mid_y = traj[:, 1].mean()
    mid_z = traj[:, 2].mean()
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    if save_png:
        plt.savefig(save_png, facecolor=fig.get_facecolor(), edgecolor="none", bbox_inches="tight")
        print(f"[Plot] Saved Matplotlib figure to: {save_png}")

    plt.show()


def plot_plotly(body_name, traj_km, earth_km, ets, unit="AU", save_html="trajectory_3d.html"):
    """
    Creates an interactive 3D trajectory visualization with time scrubbing using Plotly.
    """
    import numpy as np
    import plotly.graph_objects as go

    scale = KM_PER_AU if unit.upper() == "AU" else 1.0
    u_label = "AU" if unit.upper() == "AU" else "km"

    traj = traj_km / scale
    earth = earth_km / scale

    # Build date hover text
    hover_dates = [format_et(et) for et in ets]
    hover_traj = [
        f"<b>{body_name}</b><br>Date: {d}<br>X: {x:.3f} {u_label}<br>Y: {y:.3f} {u_label}<br>Z: {z:.3f} {u_label}"
        for d, (x, y, z) in zip(hover_dates, traj)
    ]

    num_frames = min(len(ets), 200)
    frame_indices = np.unique(np.linspace(0, len(ets) - 1, num_frames, dtype=int))
    init_idx = int(frame_indices[-1])

    fig = go.Figure()

    # 0. Sun (static)
    fig.add_trace(go.Scatter3d(
        x=[0], y=[0], z=[0],
        mode="markers",
        marker=dict(size=12, color="#ffcc00", symbol="circle", opacity=1.0),
        name="Sun",
        hoverinfo="text",
        text=["<b>Sun</b><br>Heliocentric Origin [0, 0, 0]"],
    ))

    # 1. Earth Orbit (static)
    fig.add_trace(go.Scatter3d(
        x=earth[:, 0], y=earth[:, 1], z=earth[:, 2],
        mode="lines",
        line=dict(color="#29b6f6", width=3, dash="dash"),
        name="Earth Orbit (~1 AU)",
        hoverinfo="name",
    ))

    # 2. Full Trajectory Reference Path (static)
    fig.add_trace(go.Scatter3d(
        x=traj[:, 0], y=traj[:, 1], z=traj[:, 2],
        mode="lines",
        line=dict(color="#ff3d00", width=2),
        opacity=0.35,
        name=f"{body_name} Full Path",
        hoverinfo="text",
        text=hover_traj,
    ))

    # 3. Start Point (static)
    fig.add_trace(go.Scatter3d(
        x=[traj[0, 0]], y=[traj[0, 1]], z=[traj[0, 2]],
        mode="markers",
        marker=dict(size=7, color="#00e676", symbol="circle"),
        name="Start Position",
        hoverinfo="text",
        text=[f"Start: {hover_dates[0]}"],
    ))

    # 4. End Point (static)
    fig.add_trace(go.Scatter3d(
        x=[traj[-1, 0]], y=[traj[-1, 1]], z=[traj[-1, 2]],
        mode="markers",
        marker=dict(size=8, color="#d500f9", symbol="diamond"),
        name="End Position",
        hoverinfo="text",
        text=[f"End: {hover_dates[-1]}"],
    ))

    # 5. Traversed Spacecraft Trajectory (dynamic)
    fig.add_trace(go.Scatter3d(
        x=traj[:init_idx + 1, 0], y=traj[:init_idx + 1, 1], z=traj[:init_idx + 1, 2],
        mode="lines",
        line=dict(color="#ff3d00", width=5),
        name=f"{body_name} Trajectory",
        hoverinfo="skip",
    ))

    # 6. Current Time Position Marker (dynamic)
    fig.add_trace(go.Scatter3d(
        x=[traj[init_idx, 0]], y=[traj[init_idx, 1]], z=[traj[init_idx, 2]],
        mode="markers",
        marker=dict(size=9, color="#ffea00", symbol="circle", line=dict(color="#ffffff", width=2)),
        name=f"{body_name} Current",
        hoverinfo="text",
        text=[hover_traj[init_idx]],
    ))

    # Animation frames and slider steps targeting traces 5 and 6
    frames = []
    slider_steps = []
    for step_i, idx in enumerate(frame_indices):
        frame_name = f"frame_{step_i}"
        frames.append(
            go.Frame(
                name=frame_name,
                traces=[5, 6],
                data=[
                    go.Scatter3d(
                        x=traj[:idx + 1, 0],
                        y=traj[:idx + 1, 1],
                        z=traj[:idx + 1, 2],
                    ),
                    go.Scatter3d(
                        x=[traj[idx, 0]],
                        y=[traj[idx, 1]],
                        z=[traj[idx, 2]],
                        text=[hover_traj[idx]],
                    ),
                ],
            )
        )
        slider_steps.append(
            dict(
                method="animate",
                label=hover_dates[idx],
                args=[
                    [frame_name],
                    dict(
                        mode="immediate",
                        frame=dict(duration=0, redraw=True),
                        transition=dict(duration=0),
                    ),
                ],
            )
        )

    fig.frames = frames

    # Pin 3D axis bounds so the bounding box stays fixed while scrubbing
    all_pts = np.vstack([traj, earth, np.zeros((1, 3))])

    def _padded_range(vals):
        vmin = float(vals.min())
        vmax = float(vals.max())
        span = vmax - vmin
        pad = span * 0.05 if span > 0 else 0.1
        return [vmin - pad, vmax + pad]

    xaxis_range = _padded_range(all_pts[:, 0])
    yaxis_range = _padded_range(all_pts[:, 1])
    zaxis_range = _padded_range(all_pts[:, 2])

    # Layout styling
    start_date = hover_dates[0]
    end_date = hover_dates[-1]
    fig.update_layout(
        title=f"<b>{body_name}</b> 3D Heliocentric Trajectory (ECLIPJ2000)<br><sup>{start_date} to {end_date}</sup>",
        template="plotly_dark",
        margin=dict(l=0, r=0, b=90, t=60),
        legend=dict(x=0.02, y=0.98, bgcolor="rgba(20, 24, 33, 0.75)"),
        scene=dict(
            xaxis=dict(title=f"X ({u_label})", range=xaxis_range),
            yaxis=dict(title=f"Y ({u_label})", range=yaxis_range),
            zaxis=dict(title=f"Z ({u_label})", range=zaxis_range),
            aspectmode="data",
            camera=dict(eye=dict(x=1.3, y=1.3, z=0.8)),
        ),
        updatemenus=[
            dict(
                type="buttons",
                showactive=False,
                x=0.02,
                y=0.0,
                xanchor="left",
                yanchor="top",
                direction="left",
                pad=dict(r=10, t=20),
                bgcolor="rgba(20, 24, 33, 0.85)",
                bordercolor="#37474f",
                font=dict(color="#e0e0e0"),
                buttons=[
                    dict(
                        label="▶ Play",
                        method="animate",
                        args=[
                            None,
                            dict(
                                frame=dict(duration=40, redraw=True),
                                fromcurrent=True,
                                transition=dict(duration=0),
                            ),
                        ],
                    ),
                    dict(
                        label="❚❚ Pause",
                        method="animate",
                        args=[
                            [None],
                            dict(
                                mode="immediate",
                                frame=dict(duration=0, redraw=False),
                                transition=dict(duration=0),
                            ),
                        ],
                    ),
                ],
            )
        ],
        sliders=[
            dict(
                active=len(slider_steps) - 1,
                x=0.14,
                y=0.0,
                len=0.84,
                xanchor="left",
                yanchor="top",
                pad=dict(b=10, t=20),
                currentvalue=dict(
                    prefix="Time: ",
                    visible=True,
                    xanchor="left",
                    font=dict(size=13, color="#e0e0e0"),
                ),
                steps=slider_steps,
            )
        ],
    )

    if save_html:
        fig.write_html(save_html)
        print(f"[Plot] Saved interactive Plotly visualization to: {save_html}")

    return fig


def main():
    parser = argparse.ArgumentParser(
        description="Plot 3D spacecraft / body trajectory from a SPICE SPK (.bsp / .spk) file."
    )
    parser.add_argument("spk_file", nargs="?", help="Path to SPICE SPK (.bsp / .spk) file")
    parser.add_argument("--target", help="NAIF ID or Name of target body (auto-detected if omitted)")
    parser.add_argument("--observer", default="SUN", help="Observer center (default: 'SUN')")
    parser.add_argument("--frame", default="ECLIPJ2000", help="Reference frame (default: 'ECLIPJ2000')")
    parser.add_argument("--points", type=int, default=600, help="Number of trajectory sample points (default: 600)")
    parser.add_argument("--unit", choices=["AU", "km"], default="AU", help="Distance unit (default: 'AU')")
    parser.add_argument("--renderer", choices=["both", "matplotlib", "plotly"], default="both", help="Plot renderer")
    parser.add_argument("--extra-kernel", action="append", default=[], help="Additional kernel path (.tls, .bsp)")
    parser.add_argument("--output-html", default="trajectory_3d.html", help="Output path for Plotly HTML")
    parser.add_argument("--output-png", default="trajectory_3d.png", help="Output path for Matplotlib PNG")

    args = parser.parse_args()

    # Verify dependencies before execution
    check_dependencies()

    # Prompt if no file was supplied on command line
    spk_path = args.spk_file
    if not spk_path:
        print("--- SPICE SPK 3D Trajectory Visualizer ---")
        spk_path = input("Enter path to SPK (.bsp / .spk) file: ").strip().strip('"\'')

    if not spk_path or not os.path.isfile(spk_path):
        print(f"\n[ERROR] File not found: '{spk_path}'")
        print("\nUsage example:")
        print("  python plot_spk_trajectory.py path/to/voyager1.bsp")
        print("  python plot_spk_trajectory.py path/to/kernel.bsp --target -98 --renderer both")
        print("\nYou can download free NASA SPK kernels from:")
        print("  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/")
        print("  https://ssd.jpl.nasa.gov/horizons/ (Vector Table -> SPK file)")
        sys.exit(1)

    # 1. Load SPICE Kernels
    load_spice_kernels(spk_path, extra_kernels=args.extra_kernel)

    # 2. Inspect available bodies
    bodies = inspect_spk_bodies(spk_path)
    if not bodies:
        print(f"[ERROR] No valid ephemeris coverage found in '{spk_path}'.")
        sys.exit(1)

    print("\nBodies found in SPK kernel:")
    for i, b in enumerate(bodies):
        print(f"  [{i+1}] ID: {b['id']} ({b['name']}) | Range: {format_et(b['start_et'])} -> {format_et(b['end_et'])}")

    # 3. Select Target Body
    if args.target:
        # Match by name or ID
        matched = None
        for b in bodies:
            if str(b["id"]) == str(args.target) or b["name"].upper() == str(args.target).upper():
                matched = b
                break
        if not matched:
            print(f"[WARNING] Target '{args.target}' not listed in SPK index; querying directly.")
            selected_body = {"id": args.target, "name": str(args.target), "start_et": bodies[0]["start_et"], "end_et": bodies[0]["end_et"]}
        else:
            selected_body = matched
    else:
        # Default to first body
        selected_body = bodies[0]

    print(f"\nPlotting trajectory for: {selected_body['name']} (ID: {selected_body['id']})")

    # 4. Sample Trajectory and Earth Orbit
    traj_km, ets = sample_trajectory(
        selected_body["id"],
        selected_body["start_et"],
        selected_body["end_et"],
        num_points=args.points,
        observer=args.observer,
        frame=args.frame,
    )

    earth_km = sample_earth_orbit(
        selected_body["start_et"],
        selected_body["end_et"],
        num_points=200,
        observer=args.observer,
        frame=args.frame,
    )

    # 5. Render Visualizations
    if args.renderer in ["both", "plotly"]:
        fig = plot_plotly(
            selected_body["name"],
            traj_km,
            earth_km,
            ets,
            unit=args.unit,
            save_html=args.output_html,
        )
        try:
            fig.show()
        except Exception:
            pass

    if args.renderer in ["both", "matplotlib"]:
        plot_matplotlib(
            selected_body["name"],
            traj_km,
            earth_km,
            ets,
            unit=args.unit,
            save_png=args.output_png,
        )


if __name__ == "__main__":
    main()
