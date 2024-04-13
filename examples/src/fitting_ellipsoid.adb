--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

pragma Ada_2022;
pragma Assertion_Policy (Check);

with Ada.Numerics.Elementary_Functions;

with Curve_Fit.Ellipsoid;
with Ada.Text_IO;
with Ada.Float_Text_IO;

with PLY;

procedure Fitting_Ellipsoid is
   use Ada.Numerics.Elementary_Functions;
   use all type Curve_Fit.Ellipsoid_Geometric_Parameter_Index;

   function To_Point (A, B : Float) return Curve_Fit.Ellipsoid.Vector_3D is
     (10.0 * Cos (A) * Cos (B) + Sin (A + 10.0 * B),
      8.0 * Sin (A) * Cos (B) + Sin (10.0 * A + B),
      6.0 * Sin (B) + Sin (5.0 * A + 3.0 * B));

   Params : constant Curve_Fit.Ellipsoid.Frame_Parameters :=
     (Center_X => 1.0,
      Center_Y => 2.0,
      Center_Z => 3.0,
      Roll     => 0.5,
      Pitch    => 0.7,
      Yaw      => 0.9);

   Points : constant Curve_Fit.Ellipsoid.Vector_List :=
     (for J in 1 .. 100 =>
        Curve_Fit.Ellipsoid.From_Canonical_Frame
          (To_Point (2.0 * Float (J), 101.0 * Float (J)), Params));

   function Image (F : Float) return String;

   function Image (F : Float) return String is
      Result : String (1 .. 7);
   begin
      Ada.Float_Text_IO.Put (Result, F, Aft => 4, Exp => 0);
      return Result;
   end Image;

   Result : Curve_Fit.Ellipsoid.Parameters;
   RSS    : Float;
begin
   Ada.Text_IO.Put_Line ("Fitting an ellipsoid:");

   for P of Points loop
      declare
         use Curve_Fit.Ellipsoid;
         T : constant Curve_Fit.Ellipsoid.Vector_3D :=
           Curve_Fit.Ellipsoid.To_Canonical_Frame (P, Params);
         R : constant Curve_Fit.Ellipsoid.Vector_3D :=
           Curve_Fit.Ellipsoid.From_Canonical_Frame (T, Params);
         D : constant Curve_Fit.Ellipsoid.Vector_3D := R - P;
         Sum : Float := 0.0;
      begin
         for X of D loop
            Sum := Sum + X ** 2;
         end loop;
         pragma Assert (Sum < 0.001);
      end;
   end loop;

   Curve_Fit.Ellipsoid.Ellipsoid_Fit
     (Result  => Result,
      RSS     => RSS,
      Points  => Points,
      Initial => (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, 6.0, 1.0),
      Epsilon => Float'Model_Epsilon);

   for J in Result'Range loop
      Ada.Text_IO.Put (J'Image);
      Ada.Text_IO.Put (":");
      Ada.Float_Text_IO.Put (Result (J), Aft => 5, Exp => 0);
      Ada.Text_IO.New_Line;
   end loop;

   Ada.Text_IO.Put_Line ("RSS:");
   Ada.Float_Text_IO.Put (RSS);
   Ada.Text_IO.New_Line;

   declare
      use type Curve_Fit.Ellipsoid.Vector_3D;
      Frame     : Curve_Fit.Ellipsoid.Frame_Parameters renames
        Result (Curve_Fit.Ellipsoid.Frame_Parameters'Range);
      Ellipsoid : Curve_Fit.Ellipsoid.Ellipsoid_Parameters renames
        Result (Curve_Fit.Ellipsoid.Ellipsoid_Parameters'Range);

      B : Curve_Fit.Ellipsoid.Vector_3D;
   begin
      for P of Points loop
         B := Curve_Fit.Ellipsoid.From_Canonical_Frame
           (Curve_Fit.Ellipsoid.Ellipsoid_Projection
              (Curve_Fit.Ellipsoid.To_Canonical_Frame (P, Frame),
               Ellipsoid,
               0.0001),
            Frame);

         PLY.Append_Face
           ([PLY.Point (P),
             PLY.Point (B + [0.05, 0.0, 0.0]),
             PLY.Point (B + [0.0, 0.05, 0.0])]);
         PLY.Append_Face
           ([PLY.Point (P),
             PLY.Point (B + [0.0, 0.05, 0.0]),
             PLY.Point (B + [0.0, 0.0, 0.05])]);
         PLY.Append_Face
           ([PLY.Point (P),
             PLY.Point (B + [0.0, 0.0, 0.05]),
             PLY.Point (B + [0.05, 0.0, 0.0])]);
      end loop;

      PLY.Write (Scale => 0.1);
   end;


   --  pragma Assert (Image (Result (Center_X)) = " 2.6996");
   --  pragma Assert (Image (Result (Center_Y)) = " 3.8160");
   --  pragma Assert (Image (Result (Semi_Major_Axis)) = " 6.5187");
   --  pragma Assert (Image (Result (Semi_Minor_Axis)) = " 3.0319");
end Fitting_Ellipsoid;
