--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

pragma Ada_2022;
pragma Assertion_Policy (Check);

with Conic_Fit.Ellipse;
with Ada.Text_IO;
with Ada.Float_Text_IO;

procedure Draw_Ellipse is
   use all type Conic_Fit.Ellipse_Parameter_Index;

   Points : constant Conic_Fit.Ellipse.Vector_Array :=
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

   Output : Ada.Text_IO.File_Type;
   Result : Conic_Fit.Ellipse.Parameters;
   RSS    : Float;
begin
   Ada.Text_IO.Create (Output, Name => "ellipse.gnuplot");
   Ada.Text_IO.Set_Output (Output);

   Ada.Text_IO.Put_Line ("set terminal png");
   Ada.Text_IO.Put_Line ("set output 'ellipse.png'");
   Ada.Text_IO.Put_Line ("set key left");
   Ada.Text_IO.Put_Line ("set parametric");
   Ada.Text_IO.Put_Line ("set trange [0:2*pi]");
   Ada.Text_IO.Put ("plot");

   for Step in 1 .. 5 loop
      Conic_Fit.Ellipse.Ellipse_Fit
        (Result  => Result,
         RSS     => RSS,
         Points  => Points,
         Initial => [5.0, 5.0, 0.0, 3.0, 2.0],
         Epsilon   => Float'Model_Epsilon,
         Max_Steps => Step);

      Ada.Text_IO.Put (Image (Result (Center_X)));  --  center x
      Ada.Text_IO.Put ("+");
      Ada.Text_IO.Put (Image (Result (Semi_Major_Axis)));  --  a
      Ada.Text_IO.Put ("*cos(t)*cos(");
      Ada.Text_IO.Put (Image (Result (Tilt_Angle)));  --  tilt
      Ada.Text_IO.Put (")+");
      Ada.Text_IO.Put (Image (Result (Semi_Minor_Axis)));  --  b
      Ada.Text_IO.Put ("*sin(t)*sin(");
      Ada.Text_IO.Put (Image (Result (Tilt_Angle)));  --  tilt
      Ada.Text_IO.Put ("),");

      Ada.Text_IO.Put (Image (Result (Center_Y)));  --  center y
      Ada.Text_IO.Put ("+");
      Ada.Text_IO.Put (Image (Result (Semi_Major_Axis)));  --  a
      Ada.Text_IO.Put ("*cos(t)*sin(");
      Ada.Text_IO.Put (Image (Result (Tilt_Angle)));  --  tilt
      Ada.Text_IO.Put (")+ -1*");
      Ada.Text_IO.Put (Image (Result (Semi_Minor_Axis)));  --  b
      Ada.Text_IO.Put ("*sin(t)*cos(");
      Ada.Text_IO.Put (Image (Result (Tilt_Angle)));  --  tilt
      Ada.Text_IO.Put (")");
      Ada.Text_IO.Put (" title ""Iter");
      Ada.Text_IO.Put (Step'Image);
      Ada.Text_IO.Put (""", ");
   end loop;

   Ada.Text_IO.Put_Line (" '-' using 1:2 title ""data""");

   for P of Points loop
      Ada.Text_IO.Put (Image (P (1)));
      Ada.Text_IO.Put_Line (Image (P (2)));
   end loop;

   Ada.Text_IO.Put_Line ("e");
end Draw_Ellipse;
