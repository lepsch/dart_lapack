// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void chegv_2stage(final int ITYPE, final int JOBZ, final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final int W, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final Box<int> INFO,) {
  final A = A_.dim();
  final B = B_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDA, LDB, LWORK, N;
      double               RWORK( * ), W( * );
      Complex            A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

      Complex            ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      bool               LQUERY, UPPER, WANTZ;
      String             TRANS;
      int                NEIG, LWMIN, LHTRD, LWTRD, KD, IB;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV2STAGE;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL lsame, ILAENV2STAGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CHEGST, CPOTRF, CTRMM, CTRSM, CHEEV_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      UPPER = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 );

      INFO = 0;
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1;
      } else if ( !( lsame( JOBZ, 'N' ) ) ) {
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
         KD    = ILAENV2STAGE( 1, 'CHETRD_2STAGE', JOBZ, N, -1, -1, -1 );
         IB    = ILAENV2STAGE( 2, 'CHETRD_2STAGE', JOBZ, N, KD, -1, -1 );
         LHTRD = ILAENV2STAGE( 3, 'CHETRD_2STAGE', JOBZ, N, KD, IB, -1 );
         LWTRD = ILAENV2STAGE( 4, 'CHETRD_2STAGE', JOBZ, N, KD, IB, -1 );
         LWMIN = N + LHTRD + LWTRD;
         WORK[1] = SROUNDUP_LWORK(LWMIN);

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -11;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CHEGV_2STAGE ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Form a Cholesky factorization of B.

      cpotrf(UPLO, N, B, LDB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         return;
      }

      // Transform problem to standard eigenvalue problem and solve.

      chegst(ITYPE, UPLO, N, A, LDA, B, LDB, INFO );
      cheev_2stage(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO );

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         NEIG = N;
         if (INFO > 0) NEIG = INFO - 1;
         if ( ITYPE == 1 || ITYPE == 2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N';
            } else {
               TRANS = 'C';
            }

            ctrsm('Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, B, LDB, A, LDA );

         } else if ( ITYPE == 3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**H *y

            if ( UPPER ) {
               TRANS = 'C';
            } else {
               TRANS = 'N';
            }

            ctrmm('Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, B, LDB, A, LDA );
         }
      }

      WORK[1] = SROUNDUP_LWORK(LWMIN);

      }
