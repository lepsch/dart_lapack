// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void sstevd(final int JOBZ, final int N, final int D, final int E, final Matrix<double> Z_, final int LDZ, final Array<double> WORK_, final int LWORK, final Array<int> IWORK_, final int LIWORK, final Box<int> INFO,) {
  final Z = Z_.dim();
  final WORK = WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ;
      int                INFO, LDZ, LIWORK, LWORK, N;
      int                IWORK( * );
      double               D( * ), E( * ), WORK( * ), Z( LDZ, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               LQUERY, WANTZ;
      int                ISCALE, LIWMIN, LWMIN;
      double               BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TNRM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH, SLANST, SROUNDUP_LWORK;
      // EXTERNAL lsame, SLAMCH, SLANST, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSTEDC, SSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      LQUERY = ( LWORK == -1 || LIWORK == -1 );

      INFO = 0;
      LIWMIN = 1;
      LWMIN = 1;
      if ( N > 1 && WANTZ ) {
         LWMIN = 1 + 4*N + N**2;
         LIWMIN = 3 + 5*N;
      }

      if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -6;
      }

      if ( INFO == 0 ) {
         WORK[1] = SROUNDUP_LWORK(LWMIN);
         IWORK[1] = LIWMIN;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -8;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -10;
         }
      }

      if ( INFO != 0 ) {
         xerbla('SSTEVD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         if (WANTZ) Z( 1, 1 ) = ONE;
         return;
      }

      // Get machine constants.

      SAFMIN = SLAMCH( 'Safe minimum' );
      EPS = SLAMCH( 'Precision' );
      SMLNUM = SAFMIN / EPS;
      BIGNUM = ONE / SMLNUM;
      RMIN = sqrt( SMLNUM );
      RMAX = sqrt( BIGNUM );

      // Scale matrix to allowable range, if necessary.

      ISCALE = 0;
      TNRM = SLANST( 'M', N, D, E );
      if ( TNRM > ZERO && TNRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / TNRM;
      } else if ( TNRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / TNRM;
      }
      if ( ISCALE == 1 ) {
         sscal(N, SIGMA, D, 1 );
         sscal(N-1, SIGMA, E( 1 ), 1 );
      }

      // For eigenvalues only, call SSTERF.  For eigenvalues and
      // eigenvectors, call SSTEDC.

      if ( !WANTZ ) {
         ssterf(N, D, E, INFO );
      } else {
         sstedc('I', N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if (ISCALE == 1) sscal( N, ONE / SIGMA, D, 1 );

      WORK[1] = SROUNDUP_LWORK(LWMIN);
      IWORK[1] = LIWMIN;

      }
