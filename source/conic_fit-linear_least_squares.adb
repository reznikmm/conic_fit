--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

pragma Ada_2022;

procedure Conic_Fit.Linear_Least_Squares
  (Rows : Positive;
   X    : out Element_Vector)
is
   subtype Range_1 is Positive range 1 .. Rows;
   subtype Range_2 is Positive range X'Range;

   M : Element_Matrix (Range_2, Range_2) := [others => [others => Zero]];
begin
   --  Fill upper half of M = Aᵀ * A
   for K in Range_1 loop
      for I in Range_2 loop
         for J in I .. Range_2'Last loop
            pragma Warnings
              (Off, "actuals for this call may be in wrong order");
            M (I, J) := @ + A (K, I) * A (K, J);
         end loop;
      end loop;
   end loop;

   --  Copy upper half of M to its lower half
   for I in Range_2 loop
      for J in I + 1 .. Range_2'Last loop
         M (J, I) := M (I, J);
      end loop;
   end loop;

   --  Here M is (Aᵀ * A)

   M := Inverse (M);
   --  Now M is (Aᵀ * A)⁻¹

   X := [others => Zero];

   for K in Range_1 loop
      for L in Range_2 loop
         declare
            P_L_K : Element := Zero;
         begin
            for R in Range_2 loop
               P_L_K := @ + M (R, L) * A (K, R);
            end loop;

            X (L) := @ + P_L_K * B (K);
         end;
      end loop;
   end loop;
end Conic_Fit.Linear_Least_Squares;
