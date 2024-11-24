// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void chbgv(final int JOBZ, final int UPLO, final int N, final int KA, final int KB, final Matrix<double> AB_, final int LDAB, final Matrix<double> BB_, final int LDBB, final int W, final Matrix<double> Z_, final int LDZ, final Array<double> _WORK_, final Array<double> RWORK_, final Box<int> INFO,) {
  final AB = AB_.dim();
  final BB = BB_.dim();
  final Z = Z_.dim();
  final _WORK = _WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, KA, KB, LDAB, LDBB, LDZ, N;
      double               RWORK( * ), W( * );
      Complex            AB( LDAB, * ), BB( LDBB, * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               UPPER, WANTZ;
      String             VECT;
      int                IINFO, INDE, INDWRK;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHBGST, CHBTRD, CPBSTF, CSTEQR, SSTERF, XERBLA

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      UPPER = lsame( UPLO, 'U' );

      INFO = 0;
      if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( UPPER || lsame( UPLO, 'L' ) ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( KA < 0 ) {
         INFO = -4;
      } else if ( KB < 0 || KB > KA ) {
         INFO = -5;
      } else if ( LDAB < KA+1 ) {
         INFO = -7;
      } else if ( LDBB < KB+1 ) {
         INFO = -9;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -12;
      }
      if ( INFO != 0 ) {
         xerbla('CHBGV', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Form a split Cholesky factorization of B.

      cpbstf(UPLO, N, KB, BB, LDBB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         return;
      }

      // Transform problem to standard eigenvalue problem.

      INDE = 1;
      INDWRK = INDE + N;
      chbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK, RWORK( INDWRK ), IINFO );

      // Reduce to tridiagonal form.

      if ( WANTZ ) {
         VECT = 'U';
      } else {
         VECT = 'N';
      }
      chbtrd(VECT, UPLO, N, KA, AB, LDAB, W, RWORK( INDE ), Z, LDZ, WORK, IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEQR.

      if ( !WANTZ ) {
         ssterf(N, W, RWORK( INDE ), INFO );
      } else {
         csteqr(JOBZ, N, W, RWORK( INDE ), Z, LDZ, RWORK( INDWRK ), INFO );
      }
      }
