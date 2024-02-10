import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dorbdb3(M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LWORK, M, P, Q, LDX11, LDX21;
      double             PHI(*), THETA(*);
      double             TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), X11(LDX11,*), X21(LDX21,*);
      // ..

// ====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0 ;
      double             C, S;
      int                CHILDINFO, I, ILARF, IORBDB5, LLARF, LORBDB5, LWORKMIN, LWORKOPT;
      bool               LQUERY;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARF, DLARFGP, DORBDB5, DROT, XERBLA
      // ..
      // .. External Functions ..
      //- double             DNRM2;
      // EXTERNAL DNRM2
      // ..
      // .. Intrinsic Function ..
      // INTRINSIC ATAN2, COS, MAX, SIN, SQRT

      // Test input arguments

      INFO = 0;
      LQUERY = LWORK == -1;

      if ( M < 0 ) {
         INFO = -1;
      } else if ( 2*P < M || P > M ) {
         INFO = -2;
      } else if ( Q < M-P || M-Q < M-P ) {
         INFO = -3;
      } else if ( LDX11 < max( 1, P ) ) {
         INFO = -5;
      } else if ( LDX21 < max( 1, M-P ) ) {
         INFO = -7;
      }

      // Compute workspace

      if ( INFO == 0 ) {
         ILARF = 2;
         LLARF = max( P, M-P-1, Q-1 );
         IORBDB5 = 2;
         LORBDB5 = Q-1;
         LWORKOPT = max( ILARF+LLARF-1, IORBDB5+LORBDB5-1 );
         LWORKMIN = LWORKOPT;
         WORK[1] = LWORKOPT;
         if ( LWORK < LWORKMIN && !LQUERY ) {
           INFO = -14;
         }
      }
      if ( INFO != 0 ) {
         xerbla('DORBDB3', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Reduce rows 1, ..., M-P of X11 and X21

      for (I = 1; I <= M-P; I++) {

         if ( I > 1 ) {
            drot(Q-I+1, X11(I-1,I), LDX11, X21(I,I), LDX11, C, S );
         }

         dlarfgp(Q-I+1, X21(I,I), X21(I,I+1), LDX21, TAUQ1(I) );
         S = X21(I,I);
         X21[I][I] = ONE;
         dlarf('R', P-I+1, Q-I+1, X21(I,I), LDX21, TAUQ1(I), X11(I,I), LDX11, WORK(ILARF) );
         dlarf('R', M-P-I, Q-I+1, X21(I,I), LDX21, TAUQ1(I), X21(I+1,I), LDX21, WORK(ILARF) )          C = sqrt( dnrm2( P-I+1, X11(I,I), 1 )**2 + dnrm2( M-P-I, X21(I+1,I), 1 )**2 );
         THETA[I] = ATAN2( S, C );

         dorbdb5(P-I+1, M-P-I, Q-I, X11(I,I), 1, X21(I+1,I), 1, X11(I,I+1), LDX11, X21(I+1,I+1), LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO );
         dlarfgp(P-I+1, X11(I,I), X11(I+1,I), 1, TAUP1(I) );
         if ( I < M-P ) {
            dlarfgp(M-P-I, X21(I+1,I), X21(I+2,I), 1, TAUP2(I) );
            PHI[I] = ATAN2( X21(I+1,I), X11(I,I) );
            C = COS( PHI(I) );
            S = SIN( PHI(I) );
            X21[I+1][I] = ONE;
            dlarf('L', M-P-I, Q-I, X21(I+1,I), 1, TAUP2(I), X21(I+1,I+1), LDX21, WORK(ILARF) );
         }
         X11[I][I] = ONE;
         dlarf('L', P-I+1, Q-I, X11(I,I), 1, TAUP1(I), X11(I,I+1), LDX11, WORK(ILARF) );

      }

      // Reduce the bottom-right portion of X11 to the identity matrix

      for (I = M-P + 1; I <= Q; I++) {
         dlarfgp(P-I+1, X11(I,I), X11(I+1,I), 1, TAUP1(I) );
         X11[I][I] = ONE;
         dlarf('L', P-I+1, Q-I, X11(I,I), 1, TAUP1(I), X11(I,I+1), LDX11, WORK(ILARF) );
      }

      }
