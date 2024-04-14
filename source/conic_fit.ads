--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

package Conic_Fit is
   pragma Pure;

   type Ellipse_Geometric_Parameter_Index is
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

   subtype Ellipse_Parameter_Index is Ellipse_Geometric_Parameter_Index range
     Semi_Major_Axis .. Semi_Minor_Axis;

   type Ellipsoid_Geometric_Parameter_Index is
     (Center_X,          --  x₀
      Center_Y,          --  y₀
      Center_Z,          --  z₀
      Roll,              --  α
      Pitch,             --  β
      Yaw,               --  γ
      Semi_Major_Axis,   --  a
      Semi_Middle_Axis,  --  b
      Semi_Minor_Axis);  --  c

   --  Ellipse canonical equation:
   --
   --   x̆²    y̆²    z̆²
   --  --- + --- + --- = 1
   --   a²    b²    c²
   --
   --  translated and rotated with
   --
   --   x̆ =  (x-x₀)⁢cos β⁢cos γ - (y-y₀)cos β⁢sin γ + (z-z₀)⁢sin β⁢
   --
   --   y̆ =  (x-x₀)⁢(⁢cos α⁢sin γ + sin α⁢sin β⁢cos γ))
   --      + (y-y₀)(cos α⁢cos γ - sin α⁢sin β⁢sin γ))
   --      - (z-z₀)⁢sin α⁢cos β
   --
   --   z̆ =  (x-x₀)⁢(⁢sin α⁢sin γ - cos α⁢sin β⁢cos γ))
   --      + (y-y₀)(sin α⁢cos γ + cos α⁢sin β⁢sin γ))
   --      + (z-z₀)⁢cos α⁢cos β
   --

   subtype Ellipsoid_Parameter_Index is
     Ellipsoid_Geometric_Parameter_Index range
       Semi_Major_Axis .. Semi_Minor_Axis;

   type Sphere_Geometric_Parameter_Index is
     (Center_X, --  x₀
      Center_Y, --  y₀
      Center_Z, --  z₀
      Radius);  --  r

end Conic_Fit;
