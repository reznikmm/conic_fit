--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

pragma Ada_2022;
pragma Assertion_Policy (Check);

with Curve_Fit.Ellipse;
with Ada.Text_IO;
with Ada.Float_Text_IO;

procedure Fitting_Ellipse is
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

   Ada.Text_IO.Put (Image (Result (1)));
   Ada.Text_IO.Put (Image (Result (2)));
   Ada.Text_IO.Put (Image (Result (3)));
   Ada.Text_IO.Put (Image (Result (4)));
   Ada.Text_IO.Put_Line (Image (Result (5)));
   Ada.Text_IO.Put_Line ("RSS:" & Image (RSS));

   pragma Assert (Image (Result (1)) = " 2.6996");
   pragma Assert (Image (Result (2)) = " 3.8160");
   pragma Assert (Image (Result (3)) = " 0.3596");
   pragma Assert (Image (Result (4)) = " 6.5187");
   pragma Assert (Image (Result (5)) = " 3.0319");
end Fitting_Ellipse;
