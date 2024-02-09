import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dsygvx(final int ITYPE, final String JOBZ, final String RANGE, final String UPLO, final int N,
          final Matrix<double> A, final int LDA,
          final Matrix<double> B, final int LDB, final double VL, final double VU, final int IL, final int IU, final double ABSTOL, final Box<int> M, final Array<double> W,
          final Matrix<double> Z, final int LDZ, final Array<double> WORK, final int LWORK, final Array<int> IWORK, final Array<int> IFAIL, final Box<int> INFO, ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO.value, ITYPE, IU, LDA, LDB, LDZ, LWORK, M.value, N;
      double             ABSTOL, VL, VU;
      int                IFAIL( * ), IWORK( * );
      double             A( LDA, * ), B( LDB, * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

      const              ONE = 1.0 ;
      bool               ALLEIG, INDEIG, LQUERY, UPPER, VALEIG, WANTZ;
      String             TRANS;
      int                LWKMIN, LWKOPT, NB;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      // EXTERNAL lsame, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DPOTRF, DSYEVX, DSYGST, DTRMM, DTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters.

      UPPER = lsame( UPLO, 'U' );
      WANTZ = lsame( JOBZ, 'V' );
      ALLEIG = lsame( RANGE, 'A' );
      VALEIG = lsame( RANGE, 'V' );
      INDEIG = lsame( RANGE, 'I' );
      LQUERY = ( LWORK == -1 );

      INFO.value = 0;
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO.value = -1;
      } else if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO.value = -2;
      } else if ( !( ALLEIG || VALEIG || INDEIG ) ) {
         INFO.value = -3;
      } else if ( !( UPPER || lsame( UPLO, 'L' ) ) ) {
         INFO.value = -4;
      } else if ( N < 0 ) {
         INFO.value = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO.value = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO.value = -9;
      } else {
         if ( VALEIG ) {
            if (N > 0 && VU <= VL) INFO.value = -11;
         } else if ( INDEIG ) {
            if ( IL < 1 || IL > max( 1, N ) ) {
               INFO.value = -12;
            } else if ( IU < min( N, IL ) || IU > N ) {
               INFO.value = -13;
            }
         }
      }
      if (INFO.value == 0) {
         if (LDZ < 1 || (WANTZ && LDZ < N)) {
            INFO.value = -18;
         }
      }

      if ( INFO.value == 0 ) {
         LWKMIN = max( 1, 8*N );
         NB = ilaenv( 1, 'DSYTRD', UPLO, N, -1, -1, -1 );
         LWKOPT = max( LWKMIN, ( NB + 3 )*N );
         WORK[1] = LWKOPT;

         if ( LWORK < LWKMIN && !LQUERY ) {
            INFO.value = -20;
         }
      }

      if ( INFO.value != 0 ) {
         xerbla('DSYGVX', -INFO.value );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      M.value = 0;
      if ( N == 0 ) {
         return;
      }

      // Form a Cholesky factorization of B.

      dpotrf(UPLO, N, B, LDB, INFO.value );
      if ( INFO.value != 0 ) {
         INFO.value = N + INFO.value;
         return;
      }

      // Transform problem to standard eigenvalue problem and solve.

      dsygst(ITYPE, UPLO, N, A, LDA, B, LDB, INFO.value );
      dsyevx(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M.value, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL, INFO.value );

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         if (INFO.value > 0) M.value = INFO.value - 1;
         if ( ITYPE == 1 || ITYPE == 2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N';
            } else {
               TRANS = 'T';
            }

            dtrsm('Left', UPLO, TRANS, 'Non-unit', N, M.value, ONE, B, LDB, Z, LDZ );

         } else if ( ITYPE == 3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**T*y

            if ( UPPER ) {
               TRANS = 'T';
            } else {
               TRANS = 'N';
            }

            dtrmm('Left', UPLO, TRANS, 'Non-unit', N, M.value, ONE, B, LDB, Z, LDZ );
         }
      }

      // Set WORK(1) to optimal workspace size.

      WORK[1] = LWKOPT;

      return;
      }
