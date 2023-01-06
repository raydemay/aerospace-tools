import math
from manim import *


class Title(Scene):
    def construct(self):
        # Title scene
        title = Tex(r"Orbit Determination with Gibbs' Method")
        self.play(Write(title, run_time=2))
        self.wait(2)
        self.clear()


class ConicSections(Scene):
    def construct(self):
        # Create Mobjects
        orbitCircle = Circle(radius=3.0, color=WHITE)
        orbitEllipse = Ellipse(width=5.0, height=2.0, color=GREEN_E)
        orbitParabola = ImplicitFunction(
            lambda x, y: (y**2) - 2.5 - x, color=PURE_RED, use_smoothing=True
        )
        focus = Dot(radius=0.1, color=BLUE)

        # Draw
        self.add(focus)
        self.play(Create(orbitCircle))
        self.wait()
        self.play(
            focus.animate.shift(LEFT), ReplacementTransform(orbitCircle, orbitEllipse)
        )  # ReplacementTransform replaces orbitCircle with orbitEllipse rather than just reshaping it
        self.wait()
        self.play(Transform(orbitEllipse, orbitParabola))
        self.wait()
        self.clear()


class OrbitalElements(Scene):
    def construct(self):
        # Create Mobjects
        heading = Tex(r"Orbital Elements")
        elements = Tex(r"a, e, ", r"i, $\omega$, $\Omega$, ", r"f")
        OEfromRV_axiom = Tex(r"These elements can be found with two vectors:")
        RV = MathTex(r"\vec{r}, \vec{v}")
        caption1 = Tex(r"Conic parameters", font_size=36).shift(DOWN)
        caption2 = Tex(r"Euler angles", font_size=36).shift(DOWN)
        caption3 = Tex(r"True anomaly", font_size=36).shift(DOWN)
        framebox1 = SurroundingRectangle(elements[0], buff=0.1)
        framebox2 = SurroundingRectangle(elements[1], buff=0.1)
        framebox3 = SurroundingRectangle(elements[2], buff=0.1)
        OEGroup = VGroup(
            heading, VGroup(elements, framebox1, framebox2, framebox3)
        ).arrange(DOWN)
        RVGroup = VGroup(OEfromRV_axiom, RV).arrange(DOWN)

        # Draw
        self.play(Write(heading))
        self.play(FadeIn(elements, shift=DOWN))
        self.wait()
        self.play(Create(framebox1), Create(caption1))
        self.wait()
        self.play(
            ReplacementTransform(framebox1, framebox2),
            ReplacementTransform(caption1, caption2),
        )
        self.wait()
        self.play(
            ReplacementTransform(framebox2, framebox3),
            ReplacementTransform(caption2, caption3),
        )
        self.wait()
        self.play(
            FadeOut(framebox3), FadeOut(caption3), FadeOut(elements), FadeOut(heading)
        )
        self.play(Write(RVGroup))
        self.wait(4)
        self.play(FadeOut(RVGroup))


class GibbsMethodStatement(Scene):
    def construct(self):
        # Create Mobjects
        central_body = Circle(radius=1.0, color=BLUE)
        central_body.set_fill(BLUE, opacity=1)
        space_object_path_1 = ParametricFunction(
            lambda t: np.array(
                [
                    -2 + 5 * np.cos(t),
                    3 * np.sin(t),
                    0,
                ]
            ),
            t_range=[-0.5, 0.4],
            color=GRAY,
        )
        space_object_path_2 = ParametricFunction(
            lambda t: np.array(
                [
                    -2 + 5 * np.cos(t),
                    3 * np.sin(t),
                    0,
                ]
            ),
            t_range=[0.4, 1.1],
            color=GRAY,
        )
        space_object = Dot(radius=0.08, color=RED).shift(
            space_object_path_1.get_point_from_function(-0.5)
        )
        to_object = Line(central_body.get_center(), space_object.get_center())
        from_object = Line(space_object.get_center(), central_body.get_center())
        to_object.add_updater(
            lambda z: z.become(
                Line(central_body.get_center(), space_object.get_center())
            )
        )
        from_object.add_updater(
            lambda w: w.become(
                Line(space_object.get_center(), central_body.get_center())
            )
        )
        # Draw
        self.add(central_body)
        self.wait(0.25)
        self.play(Create(space_object))
        self.wait()
        self.play(ShowPassingFlash(to_object, time_width=0.5))
        self.play(ShowPassingFlash(from_object, time_width=0.5))
        self.wait()
        self.play(MoveAlongPath(space_object, space_object_path_1))
        self.wait(0.1)
        self.play(ShowPassingFlash(to_object, time_width=0.5))
        self.play(ShowPassingFlash(from_object, time_width=0.5))
        self.wait()
        self.play(MoveAlongPath(space_object, space_object_path_2))
        self.wait(0.1)
        self.play(ShowPassingFlash(to_object, time_width=0.5))
        self.play(ShowPassingFlash(from_object, time_width=0.5))
        self.wait()
        self.clear()

        self.next_section("Method Intro")
        # Create Mobjects
        aha = Tex(r"These measurements are \emph{position} vectors")
        Earth = Dot(radius=0.2, color=BLUE)
        arrow = Tex(r"$\xrightarrow{?}$").next_to(Earth, 8 * DOWN + LEFT)
        input = Tex(r"($\vec{r_1}$, $\vec{r_2}$, $\vec{r_3}$)").next_to(arrow, LEFT)
        elements = Tex(r"(a, e, i, $\omega$, $\Omega$, f)").next_to(arrow, RIGHT)
        curve = ParametricFunction(
            lambda t: np.array(
                [
                    -2 + 5 * np.cos(t),
                    3 * np.sin(t),
                    0,
                ]
            ),
            t_range=[-0.5, 1.1],
            color=GRAY,
        )  # Parametric equation for ellipse
        r1 = Arrow(
            Earth.get_center(),
            curve.get_point_from_function(-0.5),
            stroke_width=3,
            buff=0.01,
        )
        r1text = Tex(r"$\vec{r_1}$").next_to(r1.get_end())  # r1 vector notation
        r2 = Arrow(
            Earth.get_center(),
            curve.get_point_from_function(0.4),
            stroke_width=3,
            buff=0.01,
        )
        r2text = Tex(r"$\vec{r_2}$").next_to(r2.get_end())  # r2 vector notation
        r3 = Arrow(
            Earth.get_center(),
            curve.get_point_from_function(1.1),
            stroke_width=3,
            buff=0.01,
        )
        r3text = Tex(r"$\vec{r_3}$").next_to(
            r3.get_end(), 0.9 * UP
        )  # r3 vector notation

        # Draw
        self.play(Create(aha))
        self.wait(2)
        self.play(FadeOut(aha))
        self.play(Create(curve), FadeIn(Earth))
        self.wait()
        self.play(
            Create(r1),
            Create(r2),
            Create(r3),
            Create(r1text),
            Create(r2text),
            Create(r3text),
        )
        self.wait()
        self.play(Write(input), Write(arrow), Write(elements))
        self.wait(2)
        self.clear()


class GibbsMethodSetup(Scene):
    def construct(self):
        # Create Mobjects
        axiom1 = Tex(r"Two-body orbits are \emph{planar} ")
        norm1 = MathTex(
            r"{{\hat{n}_{12}}}",
            r"\triangleq",
            r"\frac{\vec{r_1}\times\vec{r_2}}{\Vert\vec{r_1}\Vert\Vert\vec{r_2}\Vert}",
        )
        norm2 = MathTex(
            r"{{\hat{i}_h}}",
            r" = ",
            r"\frac{\vec{r_1}\times\vec{r_2}}{\Vert\vec{r_1}\Vert\Vert\vec{r_2}\Vert}",
        )
        i_h = MathTex(r"\hat{i}_h", color=GREEN).shift(3 * LEFT + DOWN)
        norm1.set_color_by_tex("12", YELLOW)
        norm2.set_color_by_tex("{i}", GREEN)
        group1 = VGroup(axiom1, norm1).arrange(
            DOWN
        )  # Group probably didn't need to be named since it isn't used
        group2 = VGroup(axiom1, norm2).arrange(DOWN)
        axiom2 = Tex(r"Two vectors form a \emph{basis} of the plane")
        linear_combo = MathTex(
            r"\vec{r_2} ",
            r" = c_1",
            r"\vec{r_1} ",
            r"+ ",
            r"c_2",
            r"\vec{r_3}",
        )
        linear_combo.set_color_by_tex("r_3", GREEN)
        linear_combo.set_color_by_tex("r_1", GREEN)
        linear_combo.set_color_by_tex("r_2", GREEN)
        group3 = VGroup(axiom2, linear_combo).arrange(DOWN)
        axiom3 = Tex(r"The motion of the velocity vector traces a \emph{hodograph} ")
        hodograph_eq = MathTex(
            r"\vec{v} = { {{\mu}} \over {h} }( e{{\hat{i}_h}}\times\hat{i}_e + {{\hat{i}_h}}\times",
            r"\frac{\vec{r}}{\Vert\vec{r}\Vert} )",
        )
        hodograph_eq.set_color_by_tex(r"{i}_h", GREEN)
        hodograph_eq.set_color_by_tex("mu", GREEN)
        hodograph_eq.set_color_by_tex(r"vec{r}", GREEN)
        group4 = VGroup(axiom3, hodograph_eq).arrange(DOWN)
        footnote = (
            Tex(r"(knowns in green)", color=GREEN).scale(0.4).to_corner(DOWN + RIGHT)
        )
        uhoh = Tex(r"How do we deal with the rest of the unknowns?")

        # Draw
        self.play(Write(axiom1))
        self.play(FadeIn(norm1))
        self.wait()
        self.play(TransformMatchingTex(Group(norm1, i_h), norm2))
        self.play(group2.animate.scale(0.8).to_corner(UP + LEFT))

        self.play(Write(axiom2))
        self.wait()
        self.play(Create(linear_combo))
        self.wait()
        self.play(group3.animate.scale(0.8).to_corner(UP + RIGHT))
        self.wait()

        self.play(Create(group4))
        self.wait()
        self.play(group4.animate.scale(0.8).to_edge(DOWN))

        self.play(FadeIn(footnote))

        self.play(Write(uhoh))
        self.wait(5)

        self.play(
            FadeOut(group2),
            FadeOut(group3),
            FadeOut(group4),
            FadeOut(footnote),
            FadeOut(uhoh),
        )


class h_and_c(Scene):
    def construct(self):
        simplified_problem = MathTex(
            r"\left(\vec{r_1}, \vec{r_2}, \vec{r_3}\right)",
            r"\rightarrow",
            r"( \vec{h}, \vec{e} )",
        )
        angular_momentum = MathTex(r"\vec{h} = \vec{r}\times\vec{v}")
        eccentricity = MathTex(
            r"{\mu}\vec{e} = \vec{v}\times\vec{h} - \frac{\mu}{\Vert\vec{r}\Vert}\vec{r}"
        )
        VGroup(angular_momentum, eccentricity).arrange(DOWN)
        simplified_problem[0].set_color(BLUE)
        simplified_problem[2].set_color(RED)
        self.play(Create(simplified_problem))
        self.wait()
        self.play(
            simplified_problem.animate.shift(2 * UP),
            Create(angular_momentum),
            Create(eccentricity),
        )
        self.wait(2)
        self.play(
            FadeOut(simplified_problem),
            FadeOut(angular_momentum),
            FadeOut(eccentricity),
        )


class AngularMomentum(Scene):
    def construct(self):
        h_setup = Tex(r"Since the position vectors are coplanar:")
        linear_combo_1 = MathTex(
            r"\vec{r_2} ",
            r" = c_1",
            r"\vec{r_1} ",
            r"+ ",
            r"c_2",
            r"\vec{r_3}",
        )
        linear_combo_2 = MathTex(
            r"\vec{e}\cdot\vec{r_2} ",
            r" = \vec{e}\cdot(" r"c_1",
            r"\vec{r_1} ",
            r"+ ",
            r"c_2",
            r"\vec{r_3}",
            r")",
        )
        linear_combo_3 = MathTex(
            r"\vec{e}\cdot\vec{r_2} ",
            r" = ",
            r"\vec{e}\cdot c_1",
            r"\vec{r_1} ",
            r"+ ",
            r"\vec{e}\cdot c_2",
            r"\vec{r_3}",
        )
        def_rDote = MathTex(r"\vec{r}\cdot\vec{e} = \frac{h^2}{\mu} - r").shift(DOWN)
        linear_combo_4 = MathTex(
            r"\frac{h^2}{\mu} - r_2 ",
            r" = ",
            r" c_1",
            r"\left(\frac{h^2}{\mu} - r_1\right) ",
            r"+ ",
            r" c_2",
            r"\left(\frac{h^2}{\mu} - r_3\right)",
        )
        self.play(Write(h_setup))
        self.play(h_setup.animate.to_corner(UP + LEFT), Create(linear_combo_1))
        self.wait()
        self.play(
            TransformMatchingTex(linear_combo_1, linear_combo_2),
        )
        self.wait()
        self.play(TransformMatchingTex(linear_combo_2, linear_combo_3))
        self.play(linear_combo_3.animate.shift(UP), Create(def_rDote))
        self.wait()
        self.play(
            FadeOut(def_rDote), TransformMatchingTex(linear_combo_3, linear_combo_4)
        )
