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

package Conic_Fit.Generic_Ellipse is
   pragma Pure;

   subtype Vector_2D is Vectors.Vector_2D;
   subtype Vector_Array is Vectors.Vector_2D_Array;

   function "+" (L, R : Vector_2D) return Vector_2D;

   function "-" (L, R : Vector_2D) return Vector_2D;

   type Parameter_Array is
     array (Ellipse_Geometric_Parameter_Index range <>) of Number;

   subtype Frame_Parameters is Parameter_Array (Center_X .. Tilt_Angle);

   function To_Canonical_Frame
     (V : Vector_2D; P : Frame_Parameters) return Vector_2D;

   function From_Canonical_Frame
     (V : Vector_2D; P : Frame_Parameters) return Vector_2D;

   subtype Ellipse_Parameters is Parameter_Array (Ellipse_Parameter_Index);

   function Ellipse_Projection
     (Point      : Vector_2D;
      Parameters : Ellipse_Parameters;
      Epsilon    : Number) return Vector_2D;

   subtype Parameters is Parameter_Array (Ellipse_Geometric_Parameter_Index);

   procedure Ellipse_Fit
     (Result    : out Parameters;
      RSS       : out Number;
      Points    : Vector_Array;
      Initial   : Parameters;
      Epsilon   : Number;
      Max_Steps : Positive := 50);

private

   function "+" (L, R : Vector_2D) return Vector_2D is
     [L (1) + R (1), L (2) + R (2)];

   function "-" (L, R : Vector_2D) return Vector_2D is
     [L (1) - R (1), L (2) - R (2)];

   function Center (P : Frame_Parameters) return Vector_2D is
     [P (Center_X), P (Center_Y)];

end Conic_Fit.Generic_Ellipse;
