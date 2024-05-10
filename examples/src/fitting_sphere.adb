--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

with Ada.Float_Text_IO;
with Ada.Numerics.Elementary_Functions;
with Ada.Text_IO;

with Conic_Fit.Sphere;

procedure Fitting_Sphere is
   use Ada.Numerics.Elementary_Functions;

   function To_Point (A, B : Float) return Conic_Fit.Sphere.Vector_3D is
     [8.0 * Cos (A) * Cos (B) + Sin (A + 10.0 * B) + 1.0,
      8.0 * Sin (A) * Cos (B) + Sin (10.0 * A + B) + 2.0,
      8.0 * Sin (B) + Sin (5.0 * A + 3.0 * B)      + 3.0];

   Points : constant Conic_Fit.Sphere.Vector_Array :=
     [for J in 1 .. 100 =>
        (To_Point (2.0 * Float (J), 101.0 * Float (J)))];

   Sphere : Conic_Fit.Sphere.Parameters;
begin
   Ada.Text_IO.Put_Line ("Fitting a sphere:");

   Conic_Fit.Sphere.Sphere_Fit (Sphere, Points);

   for J in Sphere'Range loop
      Ada.Text_IO.Put (J'Image);
      Ada.Text_IO.Put (":");
      Ada.Float_Text_IO.Put (Sphere (J), Aft => 5, Exp => 0);
      Ada.Text_IO.New_Line;
   end loop;

end Fitting_Sphere;
