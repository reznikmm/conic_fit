--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

pragma Ada_2022;

with Conic_Fit.Generic_Vectors;

generic
   type Number is private;

   One : Number;

   with function "+" (L, R : Number) return Number is <>;
   with function "-" (L, R : Number) return Number is <>;
   with function "*" (L, R : Number) return Number is <>;
   with function "/" (L, R : Number) return Number is <>;
   with function "<" (L, R : Number) return Boolean is <>;

   with function Sin (L : Number) return Number is <>;
   with function Cos (L : Number) return Number is <>;
   with function Sqrt (L : Number) return Number is <>;

   type Matrix is array (Integer range <>, Integer range <>) of Number;
   type Vector is array (Integer range <>) of Number;

   with function "*" (L, R : Matrix) return Matrix is <>;
   with function "*" (L : Matrix; R : Vector) return Vector is <>;

   with function Inverse (M : Matrix) return Matrix is <>;

   with package Vectors is new
     Conic_Fit.Generic_Vectors (Number, Vector);

package Conic_Fit.Generic_Ellipsoid is
   pragma Pure;

   subtype Vector_3D is Vectors.Vector_3D;
   subtype Vector_Array is Vectors.Vector_3D_Array;

   type Parameter_Array is
     array (Ellipsoid_Geometric_Parameter_Index range <>) of Number;

   subtype Frame_Parameters is Parameter_Array (Center_X .. Yaw);

   function To_Canonical_Frame
     (V : Vector_3D; P : Frame_Parameters) return Vector_3D;

   function From_Canonical_Frame
     (V : Vector_3D; P : Frame_Parameters) return Vector_3D;

   subtype Ellipsoid_Parameters is Parameter_Array (Ellipsoid_Parameter_Index);

   function Ellipsoid_Projection
     (Point      : Vector_3D;
      Parameters : Ellipsoid_Parameters;
      Epsilon    : Number) return Vector_3D;

   subtype Parameters is Parameter_Array (Ellipsoid_Geometric_Parameter_Index);

   procedure Ellipsoid_Fit
     (Result    : out Parameters;
      RSS       : out Number;
      Points    : Vector_Array;
      Initial   : Parameters;
      Epsilon   : Number;
      Max_Steps : Positive := 50);

end Conic_Fit.Generic_Ellipsoid;
