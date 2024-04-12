--  SPDX-FileCopyrightText: 2024 Max Reznik <reznikmm@gmail.com>
--
--  SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
----------------------------------------------------------------

package body Curve_Fit.Generic_Ellipsoid is

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

   function dP
     (V : Vector_3D;
      R : Parameters;
      P : Ellipsoid_Geometric_Parameter_Index) return Number
   is
      function F
        (CP : Vector_3D;
         R  : Parameters) return Number is
          (CP (1) ** 2 / R (Semi_Major_Axis) ** 2 +
           CP (2) ** 2 / R (Semi_Middle_Axis) ** 2 +
           CP (3) ** 2 / R (Semi_Minor_Axis) ** 2 - One);

      Frame : Frame_Parameters renames R (Frame_Parameters'Range);

      P0    : constant Number := F (To_Canonical_Frame (V, Frame), R);
      P1    : Number;
      Copy  : Parameters := R;
   begin
      Copy (P) := @ + One / ((3 * One) ** 2) ** 2;
      P1 := F (To_Canonical_Frame (V, Copy (Frame_Parameters'Range)), Copy);
      return (P0 - P1) / (R (P) - Copy (P));
   end dP;

   -----------------
   -- Ellipsoid_Fit --
   -----------------

   procedure Ellipsoid_Fit
     (Result    : out Parameters;
      RSS       : out Number;
      Points    : Vector_List;
      Initial   : Parameters;
      Epsilon   : Number;
      Max_Steps : Positive := 50)
   is
      subtype Matrix_T is Matrix (1 .. 9, Points'Range);

      Frame   : Frame_Parameters renames Result (Frame_Parameters'Range);
      Ellipsoid : Ellipsoid_Parameters renames
        Result (Ellipsoid_Parameter_Index);
   begin
      Result := Initial;
      RSS := -One;

      for Step in 1 .. Max_Steps loop
         declare
            A     : constant Number := Result (Semi_Major_Axis);
            B     : constant Number := Result (Semi_Middle_Axis);
            C     : constant Number := Result (Semi_Minor_Axis);
            Sin_A : constant Number := Sin (Result (Roll));
            Cos_A : constant Number := Cos (Result (Roll));
            Sin_B : constant Number := Sin (Result (Pitch));
            Cos_B : constant Number := Cos (Result (Pitch));
            Sin_C : constant Number := Sin (Result (Yaw));
            Cos_C : constant Number := Cos (Result (Yaw));

            function Inside_Ellipsoid (U, V, W : Number) return Boolean is
              ((U / A) ** 2 + (V / B) ** 2 + (W / C) ** 2 < One);
            --  u²/a² + v²/b² + w²/c²  < 1

            J_T : Matrix_T;
            M   : Matrix (Matrix_T'Range (1), Matrix_T'Range (1)) :=
              [others => [others => Zero]];

            Residuals : Vector (Points'Range);
         begin
            --  Build J_T and Residuals
            for K in Points'Range loop
               declare
                  P  : constant Vector_3D := Points (K);
                  CP : constant Vector_3D := To_Canonical_Frame (P, Frame);

                  UVW : constant Vector_3D :=
                    Ellipsoid_Projection (CP, Ellipsoid, Epsilon);

                  Pr : constant Vector_3D := From_Canonical_Frame (UVW, Frame);

                  CX : constant Number := Pr (1) - Result (Center_X);
                  CY : constant Number := Pr (2) - Result (Center_Y);
                  CZ : constant Number := Pr (3) - Result (Center_Z);

                  U : constant Number := UVW (1);
                  V : constant Number := UVW (2);
                  W : constant Number := UVW (3);

                  dP_dx : constant Number :=
                    2 * U * Cos_B * Cos_C / A ** 2 +
                    2 * V * (Cos_A * Sin_C + Sin_A * Sin_B * Cos_C) / B ** 2 +
                    2 * W * (Sin_A * Sin_C - Cos_A * Sin_B * Cos_C) / C ** 2;
                  --  ∂P / ∂x

                  dP_dy : constant Number :=
                    2 * U * (-Cos_B * Sin_C) / A ** 2 +
                    2 * V * (Cos_A * Cos_C - Sin_A * Sin_B * Sin_C) / B ** 2 +
                    2 * W * (Sin_A * Cos_C + Cos_A * Sin_B * Sin_C) / C ** 2;
                  --  ∂P / ∂y

                  dP_dz : constant Number :=
                    2 * U * Sin_B / A ** 2 +
                    2 * V * (-Sin_A * Cos_B) / B ** 2 +
                    2 * W * Cos_A * Cos_B / C ** 2;
                  --  ∂P / ∂z

                  XYZ   : constant Number :=
                    Sqrt (dP_dx ** 2 + dP_dy ** 2 + dP_dy ** 2);
                  --  Sqrt (dP_dx² + dP_dy² + dP_dz²)

                  dP_dx0 : constant Number := -dP_dx;
                  --  ∂P / ∂x₀

                  dP_dy0 : constant Number := -dP_dy;
                  --  ∂P / ∂y₀

                  dP_dz0 : constant Number := -dP_dz;
                  --  ∂P / ∂z₀

                  dP_droll : constant Number :=
                    2 * V / B ** 2 *
                      (CX * (Cos_A * Sin_B * Cos_C - Sin_A * Sin_C)
                       - CY * (Sin_A *  Cos_C + Cos_A * Sin_B * Sin_C)
                       - CZ * Cos_A * Cos_B) +
                    2 * W / C ** 2 *
                      (CX * (Cos_A * Sin_C + Sin_A * Sin_B * Cos_C)
                       + CY * (Cos_A * Cos_C - Sin_A * Sin_B * Sin_C)
                       - CZ * Sin_A * Cos_B);
                  --  ∂P / ∂α

                  dP_dpitch : constant Number :=
                    2 * U / A ** 2 *
                      (CX * (-Sin_B * Cos_C)
                       + CY * (Sin_B * Sin_C)
                       + CZ * Cos_B) +
                    2 * V / B ** 2 *
                      (CX * (Sin_A * Cos_B * Cos_C)
                       - CY * Sin_A * Cos_B * Sin_C
                       + CZ * Sin_A * Sin_B) +
                    2 * W / C ** 2 *
                      (-CX * Cos_A * Cos_B * Cos_C
                       + CY * Cos_A * Cos_B * Sin_C
                       - CZ * Cos_A * Sin_B);
                  --  ∂P / ∂β

                  dP_dyaw : constant Number :=
                    2 * U / A ** 2 *
                      (-CX * Cos_B * Sin_C
                       - CY * Cos_B * Cos_C) +
                    2 * V / B ** 2 *
                      (CX * (Cos_A * Cos_C - Sin_A * Sin_B * Sin_C)
                       - CY * (Cos_A * Sin_C + Sin_A * Sin_B * Cos_C)) +
                    2 * W / C ** 2 *
                      (CX * (Sin_A * Cos_C + Cos_A * Sin_B * Sin_C)
                       - CY * (Sin_A * Sin_C - Cos_A * Sin_B * Cos_C));
                  --  ∂P / ∂γ

                  dP_da : constant Number := -(2 * U ** 2 / A ** 3);
                  --  ∂P / ∂a

                  dP_db : constant Number := -(2 * V ** 2 / B ** 3);
                  --  ∂P / ∂b

                  dP_dc : constant Number := -(2 * W ** 2 / C ** 3);
                  --  ∂P / ∂c

                  RV : constant Vector_3D := P - Pr;
               begin
                  J_T (1, K) := dP_dx0 / XYZ;
                  J_T (2, K) := dP_dy0 / XYZ;
                  J_T (3, K) := dP_dz0 / XYZ;
                  J_T (4, K) := dP_droll / XYZ;
                  J_T (5, K) := dP_dpitch / XYZ;
                  J_T (6, K) := dP_dyaw / XYZ;
                  J_T (7, K) := dP_da / XYZ;
                  J_T (8, K) := dP_db / XYZ;
                  J_T (9, K) := dP_dc / XYZ;

                  for Ppp in Ellipsoid_Geometric_Parameter_Index loop
                     declare
                        Exp : Number := dP (Pr, Result, Ppp);
                        Got : constant Number := J_T
                          (1 + Ellipsoid_Geometric_Parameter_Index'Pos (Ppp),
                           K);
                     begin
                        if Exp < Zero xor Got < Zero then
                           Exp := Zero;
                        end if;
                     end;
                  end loop;

                  Residuals (K) :=
                    Sqrt (RV (1) ** 2 + RV (2) ** 2 + RV (3) ** 2);

                  if Inside_Ellipsoid (CP (1), CP (2), CP (3)) then
                     Residuals (K) := -Residuals (K);
                  end if;
               end;

               --  Fill upper half of M = J x Jᵀ
               for I in Matrix_T'Range (1) loop
                  for J in I .. Matrix_T'Last (1) loop
                     M (I, J) := @ + J_T (I, K) * J_T (J, K);
                  end loop;
               end loop;
            end loop;

            --  Copy upper half of M to lower half
            for I in Matrix_T'Range (1) loop
               for J in I + 1 .. Matrix_T'Last (1) loop
                  M (J, I) := M (I, J);
               end loop;
            end loop;

            M := Inverse (M);

            J_T := M * J_T;

            declare
               D : constant Vector (Matrix_T'Range (1)) := J_T * Residuals;
            begin
               for J in D'Range loop
                  Result (Ellipsoid_Parameter_Index'Val (J - 1)) := @ - D (J);
               end loop;
            end;

            declare
               Next_RSS : constant Number :=
                 [for R of Residuals => R ** 2]'Reduce ("+", Zero);
            begin
               exit when abs (RSS - Next_RSS) < Epsilon;
               RSS := Next_RSS;
            end;
         end;
      end loop;
   end Ellipsoid_Fit;

   ------------------------
   -- Ellipsoid_Projection --
   ------------------------

   function Ellipsoid_Projection
     (Point      : Vector_3D;
      Parameters : Ellipsoid_Parameters;
      Epsilon    : Number) return Vector_3D
   is
      CP : Vector_3D := Point;
      U  : Number renames CP (1);
      V  : Number renames CP (2);
      W  : Number renames CP (3);
      A  : Number renames Parameters (Semi_Major_Axis);
      B  : Number renames Parameters (Semi_Middle_Axis);
      C  : Number renames Parameters (Semi_Minor_Axis);

      function F (T : Number) return Number;
      function DF (T : Number) return Number;

      function F (T : Number) return Number is
         P1 : constant Number := A * U / (T + A * A);
         P2 : constant Number := B * V / (T + B * B);
         P3 : constant Number := C * W / (T + C * C);
      begin
         return P1 ** 2 + P2 ** 2 + P3 ** 2 - One;
      end F;

      function DF (T : Number) return Number is
         P1 : constant Number := 2 * (A * U) ** 2 / (T + A ** 2) ** 3;
         P2 : constant Number := 2 * (B * V) ** 2 / (T + B ** 2) ** 3;
         P3 : constant Number := 2 * (C * W) ** 2 / (T + C ** 2) ** 3;
      begin
         return -P1 - P2 - P3;
      end DF;

      Turn : array (Point'Range) of Boolean := [others => False];
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

      --  Initialize T: T₀ = max {au-a**2, bv-b**2, cw-c**2}
      declare
         T1 : constant Number := A * (U - A);
         T2 : constant Number := B * (V - B);
         T3 : constant Number := C * (W - C);
      begin
         T := (if T1 < T2 then T2 else T1);
         T := (if T < T3 then T3 else T);
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
         B ** 2 * V / (T + B ** 2),
         C ** 2 * W / (T + C ** 2)];

      --  Mirror CP back
      for J in CP'Range loop
         if Turn (J) then
            CP (J) := -CP (J);
         end if;
      end loop;

      return CP;
   end Ellipsoid_Projection;

   --------------------------
   -- From_Canonical_Frame --
   --------------------------

   function From_Canonical_Frame
     (V : Vector_3D; P : Frame_Parameters) return Vector_3D
   is
      R : constant Vector_3D :=
           [1 => V (1) * Cos (P (Pitch)) * Cos (P (Yaw))
               + V (2) * (Cos (P (Roll)) * Sin (P (Yaw))
                          + Sin (P (Roll)) * Sin (P (Pitch)) * Cos (P (Yaw)))
               + V (3) * (Sin (P (Roll)) * Sin (P (Yaw))
                          - Cos (P (Roll)) * Sin (P (Pitch)) * Cos (P (Yaw))),
            2 => -V (1) * Cos (P (Pitch)) * Sin (P (Yaw))
               + V (2) * (Cos (P (Roll)) * Cos (P (Yaw))
                          - Sin (P (Roll)) * Sin (P (Pitch)) * Sin (P (Yaw)))
               + V (3) * (Sin (P (Roll)) * Cos (P (Yaw))
                          + Cos (P (Roll)) * Sin (P (Pitch)) * Sin (P (Yaw))),
            3 => V (1) * Sin (P (Pitch))
               - V (2) * Sin (P (Roll)) * Cos (P (Pitch))
               + V (3) * Cos (P (Roll)) * Cos (P (Pitch))];
      C : constant Vector_3D := Center (P);
   begin
      return C + R;
   end From_Canonical_Frame;

   ------------------------
   -- To_Canonical_Frame --
   ------------------------

   function To_Canonical_Frame
     (V : Vector_3D; P : Frame_Parameters) return Vector_3D
   is
      R : constant Vector_3D := V - Center (P);
   begin
      return [1 => R (1) * Cos (P (Pitch)) * Cos (P (Yaw))
                 - R (2) * Cos (P (Pitch)) * Sin (P (Yaw))
                 + R (3) * Sin (P (Pitch)),
              2 => R (1) * (Cos (P (Roll)) * Sin (P (Yaw))
                          + Sin (P (Roll)) * Sin (P (Pitch)) * Cos (P (Yaw)))
                 + R (2) * (Cos (P (Roll)) * Cos (P (Yaw))
                          - Sin (P (Roll)) * Sin (P (Pitch)) * Sin (P (Yaw)))
                 - R (3) * Sin (P (Roll)) * Cos (P (Pitch)),
              3 => R (1) * (Sin (P (Roll)) * Sin (P (Yaw))
                          - Cos (P (Roll)) * Sin (P (Pitch)) * Cos (P (Yaw)))
                 + R (2) * (Sin (P (Roll)) * Cos (P (Yaw))
                          + Cos (P (Roll)) * Sin (P (Pitch)) * Sin (P (Yaw)))
                 + R (3) * Cos (P (Roll)) * Cos (P (Pitch))];
   end To_Canonical_Frame;

end Curve_Fit.Generic_Ellipsoid;