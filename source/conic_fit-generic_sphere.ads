--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

pragma Ada_2022;

generic
   type Number is private;

   Zero : Number;
   One  : Number;

   with function "+" (L, R : Number) return Number is <>;
   with function "*" (L, R : Number) return Number is <>;

   with function Sqrt (L : Number) return Number is <>;

   type Matrix is array (Integer range <>, Integer range <>) of Number;
   type Vector is array (Integer range <>) of Number;

   with function Inverse (M : Matrix) return Matrix is <>;

package Conic_Fit.Generic_Sphere is
   pragma Pure;

   subtype Vector_3D is Vector (1 .. 3);

   type Parameters is
     array (Sphere_Geometric_Parameter_Index) of Number;

   type Vector_List is array (Positive range <>) of Vector_3D;

   procedure Sphere_Fit
     (Result : out Parameters;
      Points : Vector_List);

end Conic_Fit.Generic_Sphere;
