import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dsygvd(ITYPE, JOBZ, UPLO, N, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, W, final Array<double> WORK, final int LWORK, final Array<int> IWORK, final int LIWORK, Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N;
      int                IWORK( * );
      double             A( LDA, * ), B( LDB, * ), W( * ), WORK( * );
      // ..

      double             ONE;
      const              ONE = 1.0 ;
      bool               LQUERY, UPPER, WANTZ;
      String             TRANS;
      int                LIOPT, LIWMIN, LOPT, LWMIN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL DPOTRF, DSYEVD, DSYGST, DTRMM, DTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      UPPER = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 || LIWORK == -1 );

      INFO = 0;
      if ( N <= 1 ) {
         LIWMIN = 1;
         LWMIN = 1;
      } else if ( WANTZ ) {
         LIWMIN = 3 + 5*N;
         LWMIN = 1 + 6*N + 2*N**2;
      } else {
         LIWMIN = 1;
         LWMIN = 2*N + 1;
      }
      LOPT = LWMIN;
      LIOPT = LIWMIN;
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1;
      } else if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( UPPER || lsame( UPLO, 'L' ) ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      }

      if ( INFO == 0 ) {
         WORK[1] = LOPT;
         IWORK[1] = LIOPT;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -11;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -13;
         }
      }

      if ( INFO != 0 ) {
         xerbla('DSYGVD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Form a Cholesky factorization of B.

      dpotrf(UPLO, N, B, LDB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         return;
      }

      // Transform problem to standard eigenvalue problem and solve.

      dsygst(ITYPE, UPLO, N, A, LDA, B, LDB, INFO );
      dsyevd(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO );
      LOPT = INT( max( LOPT.toDouble(), (WORK( 1 )).toDouble() ) );
      LIOPT = INT( max( LIOPT.toDouble(), (IWORK( 1 )).toDouble() ) );

      if ( WANTZ && INFO == 0 ) {

         // Backtransform eigenvectors to the original problem.

         if ( ITYPE == 1 || ITYPE == 2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N';
            } else {
               TRANS = 'T';
            }

            dtrsm('Left', UPLO, TRANS, 'Non-unit', N, N, ONE, B, LDB, A, LDA );

         } else if ( ITYPE == 3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**T*y

            if ( UPPER ) {
               TRANS = 'T';
            } else {
               TRANS = 'N';
            }

            dtrmm('Left', UPLO, TRANS, 'Non-unit', N, N, ONE, B, LDB, A, LDA );
         }
      }

      WORK[1] = LOPT;
      IWORK[1] = LIOPT;

      }
