      SUBROUTINE SGET38( RMAX, LMAX, NINFO, KNT, NIN );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, NIN;
      // ..
      // .. Array Arguments ..
      int                LMAX( 3 ), NINFO( 3 );
      REAL               RMAX( 3 );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      REAL               EPSIN;
      const              EPSIN = 5.9605e-8 ;
      int                LDT, LWORK;
      const              LDT = 20, LWORK = 2*LDT*( 10+LDT ) ;
      int                LIWORK;
      const              LIWORK = LDT*LDT ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, ISCL, ITMP, J, KMIN, M, N, NDIM;
      REAL               BIGNUM, EPS, S, SEP, SEPIN, SEPTMP, SIN, SMLNUM, STMP, TNRM, TOL, TOLIN, V, VIMIN, VMAX, VMUL, VRMIN;
      // ..
      // .. Local Arrays ..
      bool               SELECT( LDT );
      int                IPNT( LDT ), ISELEC( LDT ), IWORK( LIWORK );
      REAL               Q( LDT, LDT ), QSAV( LDT, LDT ), QTMP( LDT, LDT ), RESULT( 2 ), T( LDT, LDT ), TMP( LDT, LDT ), TSAV( LDT, LDT ), TSAV1( LDT, LDT ), TTMP( LDT, LDT ), VAL( 3 ), WI( LDT ), WITMP( LDT ), WORK( LWORK ), WR( LDT ), WRTMP( LDT );
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE;
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEHRD, SHSEQR, SHST01, SLACPY, SORGHR, SSCAL, STRSEN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL, SQRT
      // ..
      // .. Executable Statements ..

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // EPSIN = 2**(-24) = precision to which input data computed

      EPS = MAX( EPS, EPSIN );
      RMAX( 1 ) = ZERO;
      RMAX( 2 ) = ZERO;
      RMAX( 3 ) = ZERO;
      LMAX( 1 ) = 0;
      LMAX( 2 ) = 0;
      LMAX( 3 ) = 0;
      KNT = 0;
      NINFO( 1 ) = 0;
      NINFO( 2 ) = 0;
      NINFO( 3 ) = 0;

      VAL( 1 ) = SQRT( SMLNUM );
      VAL( 2 ) = ONE;
      VAL( 3 ) = SQRT( SQRT( BIGNUM ) );

      // Read input data until N=0.  Assume input eigenvalues are sorted
      // lexicographically (increasing by real part, then decreasing by
      // imaginary part)

      } // 10
      READ( NIN, FMT = * )N, NDIM;
      if (N == 0) RETURN;
      READ( NIN, FMT = * )( ISELEC( I ), I = 1, NDIM );
      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N );
      } // 20
      READ( NIN, FMT = * )SIN, SEPIN;

      TNRM = SLANGE( 'M', N, N, TMP, LDT, WORK );
      for (ISCL = 1; ISCL <= 3; ISCL++) { // 160

         // Scale input matrix

         KNT = KNT + 1;
         slacpy('F', N, N, TMP, LDT, T, LDT );
         VMUL = VAL( ISCL );
         for (I = 1; I <= N; I++) { // 30
            sscal(N, VMUL, T( 1, I ), 1 );
         } // 30
         if (TNRM == ZERO) VMUL = ONE;
         slacpy('F', N, N, T, LDT, TSAV, LDT );

         // Compute Schur form

         sgehrd(N, 1, N, T, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, INFO );
         if ( INFO != 0 ) {
            LMAX( 1 ) = KNT;
            NINFO( 1 ) = NINFO( 1 ) + 1;
            GO TO 160;
         }

         // Generate orthogonal matrix

         slacpy('L', N, N, T, LDT, Q, LDT );
         sorghr(N, 1, N, Q, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, INFO );

         // Compute Schur form

         shseqr('S', 'V', N, 1, N, T, LDT, WR, WI, Q, LDT, WORK, LWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 2 ) = KNT;
            NINFO( 2 ) = NINFO( 2 ) + 1;
            GO TO 160;
         }

         // Sort, select eigenvalues

         for (I = 1; I <= N; I++) { // 40
            IPNT( I ) = I;
            SELECT( I ) = false;
         } // 40
         scopy(N, WR, 1, WRTMP, 1 );
         scopy(N, WI, 1, WITMP, 1 );
         for (I = 1; I <= N - 1; I++) { // 60
            KMIN = I;
            VRMIN = WRTMP( I );
            VIMIN = WITMP( I );
            for (J = I + 1; J <= N; J++) { // 50
               if ( WRTMP( J ) < VRMIN ) {
                  KMIN = J;
                  VRMIN = WRTMP( J );
                  VIMIN = WITMP( J );
               }
            } // 50
            WRTMP( KMIN ) = WRTMP( I );
            WITMP( KMIN ) = WITMP( I );
            WRTMP( I ) = VRMIN;
            WITMP( I ) = VIMIN;
            ITMP = IPNT( I );
            IPNT( I ) = IPNT( KMIN );
            IPNT( KMIN ) = ITMP;
         } // 60
         for (I = 1; I <= NDIM; I++) { // 70
            SELECT( IPNT( ISELEC( I ) ) ) = true;
         } // 70

         // Compute condition numbers

         slacpy('F', N, N, Q, LDT, QSAV, LDT );
         slacpy('F', N, N, T, LDT, TSAV1, LDT );
         strsen('B', 'V', SELECT, N, T, LDT, Q, LDT, WRTMP, WITMP, M, S, SEP, WORK, LWORK, IWORK, LIWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 3 ) = KNT;
            NINFO( 3 ) = NINFO( 3 ) + 1;
            GO TO 160;
         }
         SEPTMP = SEP / VMUL;
         STMP = S;

         // Compute residuals

         shst01(N, 1, N, TSAV, LDT, T, LDT, Q, LDT, WORK, LWORK, RESULT );
         VMAX = MAX( RESULT( 1 ), RESULT( 2 ) );
         if ( VMAX > RMAX( 1 ) ) {
            RMAX( 1 ) = VMAX;
            if( NINFO( 1 ) == 0 ) LMAX( 1 ) = KNT;
         }

         // Compare condition number for eigenvalue cluster
         // taking its condition number into account

         V = MAX( TWO*REAL( N )*EPS*TNRM, SMLNUM );
         if (TNRM == ZERO) V = ONE;
         if ( V > SEPTMP ) {
            TOL = ONE;
         } else {
            TOL = V / SEPTMP;
         }
         if ( V > SEPIN ) {
            TOLIN = ONE;
         } else {
            TOLIN = V / SEPIN;
         }
         TOL = MAX( TOL, SMLNUM / EPS );
         TOLIN = MAX( TOLIN, SMLNUM / EPS );
         if ( EPS*( SIN-TOLIN ) > STMP+TOL ) {
            VMAX = ONE / EPS;
         } else if ( SIN-TOLIN > STMP+TOL ) {
            VMAX = ( SIN-TOLIN ) / ( STMP+TOL );
         } else if ( SIN+TOLIN < EPS*( STMP-TOL ) ) {
            VMAX = ONE / EPS;
         } else if ( SIN+TOLIN < STMP-TOL ) {
            VMAX = ( STMP-TOL ) / ( SIN+TOLIN );
         } else {
            VMAX = ONE;
         }
         if ( VMAX > RMAX( 2 ) ) {
            RMAX( 2 ) = VMAX;
            if( NINFO( 2 ) == 0 ) LMAX( 2 ) = KNT;
         }

         // Compare condition numbers for invariant subspace
         // taking its condition number into account

         if ( V > SEPTMP*STMP ) {
            TOL = SEPTMP;
         } else {
            TOL = V / STMP;
         }
         if ( V > SEPIN*SIN ) {
            TOLIN = SEPIN;
         } else {
            TOLIN = V / SIN;
         }
         TOL = MAX( TOL, SMLNUM / EPS );
         TOLIN = MAX( TOLIN, SMLNUM / EPS );
         if ( EPS*( SEPIN-TOLIN ) > SEPTMP+TOL ) {
            VMAX = ONE / EPS;
         } else if ( SEPIN-TOLIN > SEPTMP+TOL ) {
            VMAX = ( SEPIN-TOLIN ) / ( SEPTMP+TOL );
         } else if ( SEPIN+TOLIN < EPS*( SEPTMP-TOL ) ) {
            VMAX = ONE / EPS;
         } else if ( SEPIN+TOLIN < SEPTMP-TOL ) {
            VMAX = ( SEPTMP-TOL ) / ( SEPIN+TOLIN );
         } else {
            VMAX = ONE;
         }
         if ( VMAX > RMAX( 2 ) ) {
            RMAX( 2 ) = VMAX;
            if( NINFO( 2 ) == 0 ) LMAX( 2 ) = KNT;
         }

         // Compare condition number for eigenvalue cluster
         // without taking its condition number into account

         if ( SIN <= REAL( 2*N )*EPS && STMP <= REAL( 2*N )*EPS ) {
            VMAX = ONE;
         } else if ( EPS*SIN > STMP ) {
            VMAX = ONE / EPS;
         } else if ( SIN > STMP ) {
            VMAX = SIN / STMP;
         } else if ( SIN < EPS*STMP ) {
            VMAX = ONE / EPS;
         } else if ( SIN < STMP ) {
            VMAX = STMP / SIN;
         } else {
            VMAX = ONE;
         }
         if ( VMAX > RMAX( 3 ) ) {
            RMAX( 3 ) = VMAX;
            if( NINFO( 3 ) == 0 ) LMAX( 3 ) = KNT;
         }

         // Compare condition numbers for invariant subspace
         // without taking its condition number into account

         if ( SEPIN <= V && SEPTMP <= V ) {
            VMAX = ONE;
         } else if ( EPS*SEPIN > SEPTMP ) {
            VMAX = ONE / EPS;
         } else if ( SEPIN > SEPTMP ) {
            VMAX = SEPIN / SEPTMP;
         } else if ( SEPIN < EPS*SEPTMP ) {
            VMAX = ONE / EPS;
         } else if ( SEPIN < SEPTMP ) {
            VMAX = SEPTMP / SEPIN;
         } else {
            VMAX = ONE;
         }
         if ( VMAX > RMAX( 3 ) ) {
            RMAX( 3 ) = VMAX;
            if( NINFO( 3 ) == 0 ) LMAX( 3 ) = KNT;
         }

         // Compute eigenvalue condition number only and compare
         // Update Q

         VMAX = ZERO;
         slacpy('F', N, N, TSAV1, LDT, TTMP, LDT );
         slacpy('F', N, N, QSAV, LDT, QTMP, LDT );
         SEPTMP = -ONE;
         STMP = -ONE;
         strsen('E', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, LIWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 3 ) = KNT;
            NINFO( 3 ) = NINFO( 3 ) + 1;
            GO TO 160;
         }
         if (S != STMP) VMAX = ONE / EPS;
         IF( -ONE != SEPTMP ) VMAX = ONE / EPS;
         for (I = 1; I <= N; I++) { // 90
            for (J = 1; J <= N; J++) { // 80
               if( TTMP( I, J ) != T( I, J ) ) VMAX = ONE / EPS;
               IF( QTMP( I, J ) != Q( I, J ) ) VMAX = ONE / EPS;
            } // 80
         } // 90

         // Compute invariant subspace condition number only and compare
         // Update Q

         slacpy('F', N, N, TSAV1, LDT, TTMP, LDT );
         slacpy('F', N, N, QSAV, LDT, QTMP, LDT );
         SEPTMP = -ONE;
         STMP = -ONE;
         strsen('V', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, LIWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 3 ) = KNT;
            NINFO( 3 ) = NINFO( 3 ) + 1;
            GO TO 160;
         }
         if (-ONE != STMP) VMAX = ONE / EPS;
         IF( SEP != SEPTMP ) VMAX = ONE / EPS;
         for (I = 1; I <= N; I++) { // 110
            for (J = 1; J <= N; J++) { // 100
               if( TTMP( I, J ) != T( I, J ) ) VMAX = ONE / EPS;
               IF( QTMP( I, J ) != Q( I, J ) ) VMAX = ONE / EPS;
            } // 100
         } // 110

         // Compute eigenvalue condition number only and compare
         // Do not update Q

         slacpy('F', N, N, TSAV1, LDT, TTMP, LDT );
         slacpy('F', N, N, QSAV, LDT, QTMP, LDT );
         SEPTMP = -ONE;
         STMP = -ONE;
         strsen('E', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, LIWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 3 ) = KNT;
            NINFO( 3 ) = NINFO( 3 ) + 1;
            GO TO 160;
         }
         if (S != STMP) VMAX = ONE / EPS;
         IF( -ONE != SEPTMP ) VMAX = ONE / EPS;
         for (I = 1; I <= N; I++) { // 130
            for (J = 1; J <= N; J++) { // 120
               if( TTMP( I, J ) != T( I, J ) ) VMAX = ONE / EPS;
               IF( QTMP( I, J ) != QSAV( I, J ) ) VMAX = ONE / EPS;
            } // 120
         } // 130

         // Compute invariant subspace condition number only and compare
         // Do not update Q

         slacpy('F', N, N, TSAV1, LDT, TTMP, LDT );
         slacpy('F', N, N, QSAV, LDT, QTMP, LDT );
         SEPTMP = -ONE;
         STMP = -ONE;
         strsen('V', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, LIWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 3 ) = KNT;
            NINFO( 3 ) = NINFO( 3 ) + 1;
            GO TO 160;
         }
         if (-ONE != STMP) VMAX = ONE / EPS;
         IF( SEP != SEPTMP ) VMAX = ONE / EPS;
         for (I = 1; I <= N; I++) { // 150
            for (J = 1; J <= N; J++) { // 140
               if( TTMP( I, J ) != T( I, J ) ) VMAX = ONE / EPS;
               IF( QTMP( I, J ) != QSAV( I, J ) ) VMAX = ONE / EPS;
            } // 140
         } // 150
         if ( VMAX > RMAX( 1 ) ) {
            RMAX( 1 ) = VMAX;
            if( NINFO( 1 ) == 0 ) LMAX( 1 ) = KNT;
         }
      } // 160
      GO TO 10;

      // End of SGET38

      }
