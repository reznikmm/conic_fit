--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

generic
   type Number is private;

   type Vector is array (Integer range <>) of Number;

package Conic_Fit.Generic_Vectors is
   pragma Pure;

   subtype Vector_2D is Vector (1 .. 2);

   type Vector_2D_Array is array (Positive range <>) of Vector_2D;

   subtype Vector_3D is Vector (1 .. 3);

   type Vector_3D_Array is array (Positive range <>) of Vector_3D;

end Conic_Fit.Generic_Vectors;
