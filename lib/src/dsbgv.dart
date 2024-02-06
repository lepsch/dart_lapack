import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dsbgv(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z, LDZ, WORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, KA, KB, LDAB, LDBB, LDZ, N;
      double             AB( LDAB, * ), BB( LDBB, * ), W( * ), WORK( * ), Z( LDZ, * );
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
      // EXTERNAL DPBSTF, DSBGST, DSBTRD, DSTEQR, DSTERF, XERBLA

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
         xerbla('DSBGV', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Form a split Cholesky factorization of B.

      dpbstf(UPLO, N, KB, BB, LDBB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         return;
      }

      // Transform problem to standard eigenvalue problem.

      INDE = 1;
      INDWRK = INDE + N;
      dsbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK( INDWRK ), IINFO );

      // Reduce to tridiagonal form.

      if ( WANTZ ) {
         VECT = 'U';
      } else {
         VECT = 'N';
      }
      dsbtrd(VECT, UPLO, N, KA, AB, LDAB, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEQR.

      if ( !WANTZ ) {
         dsterf(N, W, WORK( INDE ), INFO );
      } else {
         dsteqr(JOBZ, N, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), INFO );
      }
      return;
      }
