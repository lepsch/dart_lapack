import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dggqrf(N, M, P, A, LDA, TAUA, B, LDB, TAUB, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDB, LWORK, M, N, P;
      double             A( LDA, * ), B( LDB, * ), TAUA( * ), TAUB( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LOPT, LWKOPT, NB, NB1, NB2, NB3;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQRF, DGERQF, DORMQR, XERBLA
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN

      // Test the input parameters

      INFO = 0;
      NB1 = ilaenv( 1, 'DGEQRF', ' ', N, M, -1, -1 );
      NB2 = ilaenv( 1, 'DGERQF', ' ', N, P, -1, -1 );
      NB3 = ilaenv( 1, 'DORMQR', ' ', N, M, P, -1 );
      NB = max( NB1, NB2, NB3 );
      LWKOPT = max( 1, max( N, M, P )*NB );
      WORK[1] = LWKOPT;
      LQUERY = ( LWORK == -1 );
      if ( N < 0 ) {
         INFO = -1;
      } else if ( M < 0 ) {
         INFO = -2;
      } else if ( P < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      } else if ( LWORK < max( 1, N, M, P ) && !LQUERY ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('DGGQRF', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // QR factorization of N-by-M matrix A: A = Q*R

      dgeqrf(N, M, A, LDA, TAUA, WORK, LWORK, INFO );
      LOPT = INT( WORK( 1 ) );

      // Update B := Q**T*B.

      dormqr('Left', 'Transpose', N, P, min( N, M ), A, LDA, TAUA, B, LDB, WORK, LWORK, INFO );
      LOPT = max( LOPT, INT( WORK( 1 ) ) );

      // RQ factorization of N-by-P matrix B: B = T*Z.

      dgerqf(N, P, B, LDB, TAUB, WORK, LWORK, INFO );

      WORK[1] = max( LOPT, INT( WORK( 1 ) ) );

      return;
      }
