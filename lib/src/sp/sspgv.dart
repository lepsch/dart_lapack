// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void sspgv(final int ITYPE, final int JOBZ, final int UPLO, final int N, final int AP, final int BP, final int W, final Matrix<double> Z_, final int LDZ, final Array<double> _WORK_, final Box<int> INFO,) {
  final Z = Z_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDZ, N;
      double               AP( * ), BP( * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               UPPER, WANTZ;
      String             TRANS;
      int                J, NEIG;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SPPTRF, SSPEV, SSPGST, STPMV, STPSV, XERBLA

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      UPPER = lsame( UPLO, 'U' );

      INFO = 0;
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1;
      } else if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( UPPER || lsame( UPLO, 'L' ) ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('SSPGV ', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Form a Cholesky factorization of B.

      spptrf(UPLO, N, BP, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         return;
      }

      // Transform problem to standard eigenvalue problem and solve.

      sspgst(ITYPE, UPLO, N, AP, BP, INFO );
      sspev(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO );

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         NEIG = N;
         if (INFO > 0) NEIG = INFO - 1;
         if ( ITYPE == 1 || ITYPE == 2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N';
            } else {
               TRANS = 'T';
            }

            for (J = 1; J <= NEIG; J++) { // 10
               stpsv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
            } // 10

         } else if ( ITYPE == 3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**T*y

            if ( UPPER ) {
               TRANS = 'T';
            } else {
               TRANS = 'N';
            }

            for (J = 1; J <= NEIG; J++) { // 20
               stpmv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
            } // 20
         }
      }
      }
