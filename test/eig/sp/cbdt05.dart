// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cbdt05(final int M, final int N, final Matrix<double> A_, final int LDA, final int S, final int NS, final Matrix<double> U_, final int LDU, final Matrix<double> VT_, final int LDVT, final Array<double> _WORK_, final int RESID,) {
  final A = A_.dim();
  final U = U_.dim();
  final VT = VT_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDU, LDVT, M, N, NS;
      double               RESID;
      double               S( * );
      Complex            A( LDA, * ), U( * ), VT( LDVT, * ), WORK( * );
      // ..

// ======================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                I, J;
      double               ANORM, EPS;
      double               DUM( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ISAMAX;
      //- REAL               SASUM, SCASUM, SLAMCH, CLANGE;
      // EXTERNAL lsame, ISAMAX, SASUM, SCASUM, SLAMCH, CLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, MAX, MIN

      // Quick return if possible.

      RESID = ZERO;
      if( min( M, N ) <= 0 || NS <= 0 ) return;

      EPS = SLAMCH( 'Precision' );
      ANORM = CLANGE( 'M', M, N, A, LDA, DUM );

      // Compute U' * A * V.

      cgemm('N', 'C', M, NS, N, CONE, A, LDA, VT, LDVT, CZERO, WORK( 1+NS*NS ), M )       CALL CGEMM( 'C', 'N', NS, NS, M, -CONE, U, LDU, WORK( 1+NS*NS ), M, CZERO, WORK, NS );

      // norm(S - U' * B * V)

      J = 0;
      for (I = 1; I <= NS; I++) { // 10
         WORK[J+I] = WORK( J+I ) + CMPLX( S( I ), ZERO );
         RESID = max( RESID, SCASUM( NS, WORK( J+1 ), 1 ) );
         J = J + NS;
      } // 10

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         if ( ANORM >= RESID ) {
            RESID = ( RESID / ANORM ) / ( REAL( N )*EPS );
         } else {
            if ( ANORM < ONE ) {
               RESID = ( min( RESID, double( N )*ANORM ) / ANORM ) / ( REAL( N )*EPS );
            } else {
               RESID = min( RESID / ANORM, REAL( N ) ) / ( REAL( N )*EPS );
            }
         }
      }

      }
