--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

pragma Ada_2022;

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

package Curve_Fit.Generic_Ellipsoid is
   pragma Pure;

   type Vector_3D is array (1 .. 3) of Number;

   function "+" (L, R : Vector_3D) return Vector_3D;

   function "-" (L, R : Vector_3D) return Vector_3D;

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
   type Vector_List is array (Positive range <>) of Vector_3D;

   procedure Ellipsoid_Fit
     (Result    : out Parameters;
      RSS       : out Number;
      Points    : Vector_List;
      Initial   : Parameters;
      Epsilon   : Number;
      Max_Steps : Positive := 50);

private

   function "+" (L, R : Vector_3D) return Vector_3D is
     [L (1) + R (1), L (2) + R (2), L (3) + R (3)];

   function "-" (L, R : Vector_3D) return Vector_3D is
     [L (1) - R (1), L (2) - R (2), L (3) - R (3)];

   function Center (P : Frame_Parameters) return Vector_3D is
     [P (Center_X), P (Center_Y), P (Center_Z)];

end Curve_Fit.Generic_Ellipsoid;
