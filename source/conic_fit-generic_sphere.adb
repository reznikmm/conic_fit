--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

with Conic_Fit.Linear_Least_Squares;

package body Conic_Fit.Generic_Sphere is

   subtype Small_Int is Positive range 2 .. 2;

   function "**" (L : Number; Ignore : Small_Int) return Number is (L * L);
   function "*" (Ignore : Small_Int; R : Number) return Number is (R + R);

   ----------------
   -- Sphere_Fit --
   ----------------

   procedure Sphere_Fit
     (Result : out Parameters;
      Points : Vector_Array)
   is

      function A (J, K : Positive) return Number is
        (case K is
            when 1 .. 3 => 2 * Points (J) (K),
            when others => One);

      function B (J : Positive) return Number is
         (Points (J) (1) ** 2 + Points (J) (2) ** 2 + Points (J) (3) ** 2);

      procedure Solve is new Conic_Fit.Linear_Least_Squares
        (Element => Number,
         Element_Vector => Vector,
         Element_Matrix => Matrix,
         Zero           => Zero,
         A              => A,
         B              => B);

      X : Vector (1 .. 4);
   begin
      Solve (Points'Last, X);

      Result :=
        [Center_X => X (1),
         Center_Y => X (2),
         Center_Z => X (3),
         Radius   => Sqrt (X (4) + X (1) ** 2 + X (2) ** 2 + X (3) ** 2)];
   end Sphere_Fit;

end Conic_Fit.Generic_Sphere;
