import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgemqrt(final int SIDE, final int TRANS, final int M, final int N, final int K, final int NB, final Matrix<double> V_, final int LDV, final Matrix<double> T_, final int LDT, final Matrix<double> C_, final int LDC, final Array<double> _WORK_, final Box<int> INFO,) {
  final V = V_.dim();
  final T = T_.dim();
  final C = C_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String    SIDE, TRANS;
      int       INFO, K, LDV, LDC, M, N, NB, LDT;
      double           V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * );
      // ..

// =====================================================================

      bool               LEFT, RIGHT, TRAN, NOTRAN;
      int                I, IB, LDWORK, KF, Q;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, DLARFB
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // .. Test the input arguments ..

      INFO   = 0;
      LEFT   = lsame( SIDE,  'L' );
      RIGHT  = lsame( SIDE,  'R' );
      TRAN   = lsame( TRANS, 'T' );
      NOTRAN = lsame( TRANS, 'N' );

      if ( LEFT ) {
         LDWORK = max( 1, N );
         Q = M;
      } else if ( RIGHT ) {
         LDWORK = max( 1, M );
         Q = N;
      }
      if ( !LEFT && !RIGHT ) {
         INFO = -1;
      } else if ( !TRAN && !NOTRAN ) {
         INFO = -2;
      } else if ( M < 0 ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( K < 0 || K > Q ) {
         INFO = -5;
      } else if ( NB < 1 || (NB > K && K > 0)) {
         INFO = -6;
      } else if ( LDV < max( 1, Q ) ) {
         INFO = -8;
      } else if ( LDT < NB ) {
         INFO = -10;
      } else if ( LDC < max( 1, M ) ) {
         INFO = -12;
      }

      if ( INFO != 0 ) {
         xerbla('DGEMQRT', -INFO );
         return;
      }

      // .. Quick return if possible ..

      if (M == 0 || N == 0 || K == 0) return;

      if ( LEFT && TRAN ) {

         for (I = 1; NB < 0 ? I >= K : I <= K; I += NB) {
            IB = min( NB, K-I+1 );
            dlarfb('L', 'T', 'F', 'C', M-I+1, N, IB, V( I, I ), LDV, T( 1, I ), LDT, C( I, 1 ), LDC, WORK, LDWORK );
         }

      } else if ( RIGHT && NOTRAN ) {

         for (I = 1; NB < 0 ? I >= K : I <= K; I += NB) {
            IB = min( NB, K-I+1 );
            dlarfb('R', 'N', 'F', 'C', M, N-I+1, IB, V( I, I ), LDV, T( 1, I ), LDT, C( 1, I ), LDC, WORK, LDWORK );
         }

      } else if ( LEFT && NOTRAN ) {

         KF = ((K-1)/NB)*NB+1;
         for (I = KF; -NB < 0 ? I >= 1 : I <= 1; I += -NB) {
            IB = min( NB, K-I+1 );
            dlarfb('L', 'N', 'F', 'C', M-I+1, N, IB, V( I, I ), LDV, T( 1, I ), LDT, C( I, 1 ), LDC, WORK, LDWORK );
         }

      } else if ( RIGHT && TRAN ) {

         KF = ((K-1)/NB)*NB+1;
         for (I = KF; -NB < 0 ? I >= 1 : I <= 1; I += -NB) {
            IB = min( NB, K-I+1 );
            dlarfb('R', 'T', 'F', 'C', M, N-I+1, IB, V( I, I ), LDV, T( 1, I ), LDT, C( 1, I ), LDC, WORK, LDWORK );
         }

      }

      }
