import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlatsqr(M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N, MB, NB, LDT, LWORK;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), WORK( * ), T( LDT, * );
      // ..

// =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, II, KK, CTR, MINMN, LWMIN;
      // ..
      // .. EXTERNAL FUNCTIONS ..
      //- bool               lsame;
      // EXTERNAL lsame
      // .. EXTERNAL SUBROUTINES ..
      // EXTERNAL DGEQRT, DTPQRT, XERBLA
      // .. INTRINSIC FUNCTIONS ..
      // INTRINSIC MAX, MIN, MOD
      // ..
      // .. EXECUTABLE STATEMENTS ..

      // TEST THE INPUT ARGUMENTS

      INFO = 0;

      LQUERY = ( LWORK == -1 );

      MINMN = min( M, N );
      if ( MINMN == 0 ) {
        LWMIN = 1;
      } else {
        LWMIN = N*NB;
      }

      if ( M < 0 ) {
        INFO = -1;
      } else if ( N < 0 || M < N ) {
        INFO = -2;
      } else if ( MB < 1 ) {
        INFO = -3;
      } else if ( NB < 1 || ( NB > N && N > 0 ) ) {
        INFO = -4;
      } else if ( LDA < max( 1, M ) ) {
        INFO = -6;
      } else if ( LDT < NB ) {
        INFO = -8;
      } else if ( LWORK < LWMIN && ( !LQUERY) ) {
        INFO = -10;
      }

      if ( INFO == 0 ) {
        WORK[1] = LWMIN;
      }

      if ( INFO != 0 ) {
        xerbla('DLATSQR', -INFO );
        return;
      } else if ( LQUERY ) {
        return;
      }

      // Quick return if possible

      if ( MINMN == 0 ) {
        return;
      }

      // The QR Decomposition

      if ( (MB <= N) || (MB >= M) ) {
        dgeqrt(M, N, NB, A, LDA, T, LDT, WORK, INFO );
        return;
      }

      KK = (M-N % MB-N);
      II = M-KK+1;

      // Compute the QR factorization of the first block A(1:MB,1:N)

      dgeqrt(MB, N, NB, A(1,1), LDA, T, LDT, WORK, INFO );

      CTR = 1;
      for (I = MB+1; (MB-N) < 0 ? I >= II-MB+N : I <= II-MB+N; I += (MB-N)) {

        // Compute the QR factorization of the current block A(I:I+MB-N,1:N)

        dtpqrt(MB-N, N, 0, NB, A(1,1), LDA, A( I, 1 ), LDA, T(1, CTR * N + 1), LDT, WORK, INFO );
        CTR = CTR + 1;
      }

      // Compute the QR factorization of the last block A(II:M,1:N)

      if ( II <= M ) {
        dtpqrt(KK, N, 0, NB, A(1,1), LDA, A( II, 1 ), LDA, T(1, CTR * N + 1), LDT, WORK, INFO );
      }

      WORK[1] = LWMIN;
      return;
      }
