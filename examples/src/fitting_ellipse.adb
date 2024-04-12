--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

pragma Ada_2022;
pragma Assertion_Policy (Check);

with Curve_Fit.Ellipse;
with Ada.Text_IO;
with Ada.Float_Text_IO;

procedure Fitting_Ellipse is
   use all type Curve_Fit.Ellipse_Geometric_Parameter_Index;

   Points : constant Curve_Fit.Ellipse.Vector_List :=
     [[1.0, 7.0],
      [2.0, 6.0],
      [5.0, 8.0],
      [7.0, 7.0],
      [9.0, 5.0],
      [3.0, 7.0],
      [6.0, 2.0],
      [8.0, 4.0]];

   function Image (F : Float) return String;

   function Image (F : Float) return String is
      Result : String (1 .. 7);
   begin
      Ada.Float_Text_IO.Put (Result, F, Aft => 4, Exp => 0);
      return Result;
   end Image;

   Result : Curve_Fit.Ellipse.Parameters;
   RSS    : Float;
begin
   Ada.Text_IO.Put ("Fitting an ellipse:");

   Curve_Fit.Ellipse.Ellipse_Fit
     (Result  => Result,
      RSS     => RSS,
      Points  => Points,
      Initial => [5.0, 5.0, 0.0, 3.0, 2.0],
      Epsilon => Float'Model_Epsilon);

   Ada.Text_IO.Put (Image (Result (Center_X)));
   Ada.Text_IO.Put (Image (Result (Center_Y)));
   Ada.Text_IO.Put (Image (Result (Tilt_Angle)));
   Ada.Text_IO.Put (Image (Result (Semi_Major_Axis)));
   Ada.Text_IO.Put_Line (Image (Result (Semi_Minor_Axis)));
   Ada.Text_IO.Put_Line ("RSS:" & Image (RSS));

   pragma Assert (Image (Result (Center_X)) = " 2.6996");
   pragma Assert (Image (Result (Center_Y)) = " 3.8160");
   pragma Assert (Image (Result (Tilt_Angle)) = " 0.3596");
   pragma Assert (Image (Result (Semi_Major_Axis)) = " 6.5187");
   pragma Assert (Image (Result (Semi_Minor_Axis)) = " 3.0319");
end Fitting_Ellipse;
