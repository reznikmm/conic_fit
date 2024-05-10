--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

pragma Ada_2022;

with Conic_Fit.Generic_Vectors;

generic
   type Number is private;

   Zero : Number;  --  0.0
   Half : Number;  --  0.5

   with function "+" (L, R : Number) return Number is <>;
   with function "*" (L, R : Number) return Number is <>;

   with function Sqrt (L : Number) return Number is <>;

   type Matrix is array (Integer range <>, Integer range <>) of Number;
   type Vector is array (Integer range <>) of Number;

   with function Inverse (M : Matrix) return Matrix is <>;

   with package Vectors is new
     Conic_Fit.Generic_Vectors (Number, Vector);

package Conic_Fit.Generic_Sphere is
   pragma Pure;

   type Parameters is
     array (Sphere_Geometric_Parameter_Index) of Number;

   subtype Vector_3D is Vectors.Vector_3D;
   subtype Vector_Array is Vectors.Vector_3D_Array;

   procedure Sphere_Fit
     (Result : out Parameters;
      Points : Vector_Array);

end Conic_Fit.Generic_Sphere;
