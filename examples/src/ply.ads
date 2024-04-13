--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

with Ada.Numerics.Real_Arrays;

package PLY is
   subtype Point is Ada.Numerics.Real_Arrays.Real_Vector (1 .. 3);
   type Face is array (1 .. 3) of Point;

   procedure Append_Face (F : Face);
   procedure Write
     (Scale : Float := 1.0;
      Name  : String := "/tmp/aaa.ply");

end PLY;
