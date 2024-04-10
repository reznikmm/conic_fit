--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

package Curve_Fit is
   pragma Pure;

   type Geometric_Parameter_Index is
     (Center_X,          --  x₀
      Center_Y,          --  y₀
      Tilt_Angle,        --  θ
      Semi_Major_Axis,   --  a
      Semi_Minor_Axis);  --  b

   --  Ellipse canonical equation:
   --
   --   x̆²    y̆²
   --  --- + --- = 1
   --   a²    b²
   --
   --  translated and rotated with
   --   x̆ =  (x-x₀)cos θ + (y-y₀)sin θ
   --   y̆ = -(x-x₀)sin θ + (y-y₀)cos θ
   --

   subtype Frame_Parameter_Index is Geometric_Parameter_Index range
     Center_X .. Tilt_Angle;

   subtype Ellipse_Parameter_Index is Geometric_Parameter_Index range
     Semi_Major_Axis .. Semi_Minor_Axis;

end Curve_Fit;
