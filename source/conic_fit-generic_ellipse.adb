--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

pragma Ada_2022;

package body Conic_Fit.Generic_Ellipse is

   function Zero return Number is (One - One);

   function "-" (L : Number) return Number is (Zero - L);

   function "abs" (L : Number) return Number is
      (if L < Zero then -L else L);

   subtype Small_Int is Positive range 2 .. 3;

   function "**" (L : Number; R : Small_Int) return Number;
   function "*" (L : Small_Int; R : Number) return Number;

   ---------
   -- "*" --
   ---------

   function "*" (L : Small_Int; R : Number) return Number is
      R_2 : constant Number := R + R;
   begin
      return (if L = 2 then R_2 else R_2 + R);
   end "*";

   ----------
   -- "**" --
   ----------

   function "**" (L : Number; R : Small_Int) return Number is
      L_2 : constant Number := L * L;
   begin
      return (if R = 2 then L_2 else L_2 * L);
   end "**";

   -----------------
   -- Ellipse_Fit --
   -----------------

   procedure Ellipse_Fit
     (Result    : out Parameters;
      RSS       : out Number;
      Points    : Vector_Array;
      Initial   : Parameters;
      Epsilon   : Number;
      Max_Steps : Positive := 50)
   is
      subtype Parameter_Index is Positive range 1 .. 5;
      --  Ellipse parameter index as an Integer range

      Frame   : Frame_Parameters renames Result (Frame_Parameters'Range);
      Ellipse : Ellipse_Parameters renames Result (Ellipse_Parameter_Index);
   begin
      Result := Initial;
      RSS := -One;

      for Step in 1 .. Max_Steps loop
         declare
            A     : constant Number := Result (Semi_Major_Axis);
            B     : constant Number := Result (Semi_Minor_Axis);
            T     : constant Number := Result (Tilt_Angle);
            Sin_T : constant Number := Sin (T);
            Cos_T : constant Number := Cos (T);

            function Inside_Ellipse (U, V : Number) return Boolean is
              ((B * U) ** 2 + (A * V) ** 2 < (A * B) ** 2);
            --  u²/a² + v²/b² < 1

            M : Matrix (Parameter_Index, Parameter_Index) :=
              [others => [others => Zero]];

            J_T_R : Vector (Parameter_Index) := [others => Zero];
            --  Jᵀ x Residuals

            Next_RSS : Number := Zero;

         begin
            --  Build M = Jᵀ x J and J_T_R = Jᵀ x Residuals
            for K in Points'Range loop
               declare
                  P  : constant Vector_2D := Points (K);
                  CP : constant Vector_2D := To_Canonical_Frame (P, Frame);

                  UV : constant Vector_2D :=
                    Ellipse_Projection (CP, Ellipse, Epsilon);

                  Pr : constant Vector_2D := From_Canonical_Frame (UV, Frame);

                  U : constant Number := UV (1);
                  V : constant Number := UV (2);

                  UB2 : constant Number := U * B ** 2;
                  VA2 : constant Number := V * A ** 2;

                  --  dP_dx : constant Number :=
                  --    2 * (Cos_T * U / A ** 2 - Sin_T * V / B ** 2);
                  --  ∂P / ∂x

                  --  dP_dy : constant Number :=
                  --    2 * (Sin_T * U / A ** 2 + Cos_T * V / B ** 2);
                  --  ∂P / ∂y

                  XY : constant Number := Sqrt (UB2 ** 2 + VA2 ** 2);
                  --  a²b² * Sqrt (dP_dx² + dP_dy²)

                  dP_dx0 : constant Number := 2 * (-Cos_T * UB2 + Sin_T * VA2);
                  --  a²b² * ∂P / ∂x₀

                  dP_dy0 : constant Number := 2 * (-Sin_T * UB2 - Cos_T * VA2);
                  --  a²b² * ∂P / ∂y₀

                  dP_dt : constant Number := 2 * (V * UB2 - U * VA2);
                  --  a²b² * ∂P / ∂θ

                  dP_da : constant Number := -(2 * U ** 2 * B ** 2 / A);
                  --  a²b² * ∂P / ∂a

                  dP_db : constant Number := -(2 * V ** 2 * A ** 2 / B);
                  --  a²b² * ∂P / ∂b

                  RV : constant Vector_2D := P - Pr;

                  Residual : Number := Sqrt (RV (1) ** 2 + RV (2) ** 2);

                  J_T : Vector (Parameter_Index);  --  K row of Jᵀ

               begin
                  J_T (1) := dP_dx0 / XY;
                  J_T (2) := dP_dy0 / XY;
                  J_T (3) := dP_dt / XY;
                  J_T (4) := dP_da / XY;
                  J_T (5) := dP_db / XY;

                  if Inside_Ellipse (CP (1), CP (2)) then
                     Residual := -Residual;
                  end if;

                  Next_RSS := @ + Residual ** 2;

                  for J in J_T_R'Range loop
                     J_T_R (J) := @ + J_T (J) * Residual;
                  end loop;

                  --  Fill upper half of M = J x Jᵀ
                  for I in Parameter_Index loop
                     for J in I .. Parameter_Index'Last loop
                        M (I, J) := @ + J_T (I) * J_T (J);
                     end loop;
                  end loop;
               end;
            end loop;

            --  Copy upper half of M to lower half
            for I in Parameter_Index loop
               for J in I + 1 .. Parameter_Index'Last loop
                  M (J, I) := M (I, J);
               end loop;
            end loop;

            M := Inverse (M);

            declare
               D : constant Vector (Parameter_Index) := M * J_T_R;
               --  D = M⁻¹ x Jᵀ x Residuals
            begin
               for J in D'Range loop
                  Result (Ellipse_Parameter_Index'Val (J - 1)) := @ - D (J);
               end loop;
            end;

            exit when abs (RSS - Next_RSS) < Epsilon;
            RSS := Next_RSS;
         end;
      end loop;
   end Ellipse_Fit;

   ------------------------
   -- Ellipse_Projection --
   ------------------------

   function Ellipse_Projection
     (Point      : Vector_2D;
      Parameters : Ellipse_Parameters;
      Epsilon    : Number) return Vector_2D
   is
      CP : Vector_2D := Point;
      U  : Number renames CP (1);
      V  : Number renames CP (2);
      A  : Number renames Parameters (Semi_Major_Axis);
      B  : Number renames Parameters (Semi_Minor_Axis);

      function F (T : Number) return Number;
      function DF (T : Number) return Number;

      function F (T : Number) return Number is
         P1 : constant Number := A * U / (T + A * A);
         P2 : constant Number := B * V / (T + B * B);
      begin
         return P1 ** 2 + P2 ** 2 - One;
      end F;

      function DF (T : Number) return Number is
         P1 : constant Number := 2 * (A * U) ** 2 / (T + A ** 2) ** 3;
         P2 : constant Number := 2 * (B * V) ** 2 / (T + B ** 2) ** 3;
      begin
         return -P1 - P2;
      end DF;

      Turn : array (Point'Range) of Boolean := [False, False];
      T    : Number;
      S    : Number;
   begin
      --  Mirror CP to the first quadrant
      for J in CP'Range loop
         if CP (J) < Zero then
            CP (J) := -CP (J);
            Turn (J) := True;
         end if;
      end loop;

      --  Initialize T: T₀ = max {au-a**2, bv-b**2}
      declare
         T1 : constant Number := A * (U - A);
         T2 : constant Number := B * (V - B);
      begin
         T := (if T1 < T2 then T2 else T1);
      end;

      --  Find T with Newton method
      for J in 1 .. 10 loop
         S := F (T) / DF (T);
         T := T - S;

         exit when abs S < Epsilon;
      end loop;

      --  Calculate the point from T
      CP :=
        [A ** 2 * U / (T + A ** 2),
         B ** 2 * V / (T + B ** 2)];

      --  Mirror CP back
      for J in CP'Range loop
         if Turn (J) then
            CP (J) := -CP (J);
         end if;
      end loop;

      return CP;
   end Ellipse_Projection;

   --------------------------
   -- From_Canonical_Frame --
   --------------------------

   function From_Canonical_Frame
     (V : Vector_2D; P : Frame_Parameters) return Vector_2D
   is
      R : constant Vector_2D :=
        [1 => V (1) * Cos (P (Tilt_Angle)) -
              V (2) * Sin (P (Tilt_Angle)),
         2 => V (1) * Sin (P (Tilt_Angle)) +
              V (2) * Cos (P (Tilt_Angle))];
      C : constant Vector_2D := Center (P);
   begin
      return C + R;
   end From_Canonical_Frame;

   ------------------------
   -- To_Canonical_Frame --
   ------------------------

   function To_Canonical_Frame
     (V : Vector_2D; P : Frame_Parameters) return Vector_2D
   is
      R : constant Vector_2D := V - Center (P);
   begin
      return [1 => R (1) * Cos (P (Tilt_Angle)) +
                   R (2) * Sin (P (Tilt_Angle)),
              2 => R (2) * Cos (P (Tilt_Angle)) -
                   R (1) * Sin (P (Tilt_Angle))];
   end To_Canonical_Frame;

end Conic_Fit.Generic_Ellipse;
