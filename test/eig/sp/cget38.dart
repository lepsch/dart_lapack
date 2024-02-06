      void cget38(RMAX, LMAX, NINFO, KNT, NIN ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                KNT, NIN;
      int                LMAX( 3 ), NINFO( 3 );
      double               RMAX( 3 );
      // ..

      int                LDT, LWORK;
      const              LDT = 20, LWORK = 2*LDT*( 10+LDT ) ;
      double               ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      double               EPSIN;
      const              EPSIN = 5.9605e-8 ;
      Complex            CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      int                I, INFO, ISCL, ISRT, ITMP, J, KMIN, M, N, NDIM;
      double               BIGNUM, EPS, S, SEP, SEPIN, SEPTMP, SIN, SMLNUM, STMP, TNRM, TOL, TOLIN, V, VMAX, VMIN, VMUL;
      bool               SELECT( LDT );
      int                IPNT( LDT ), ISELEC( LDT );
      double               RESULT( 2 ), RWORK( LDT ), VAL( 3 ), WSRT( LDT )       Complex            Q( LDT, LDT ), QSAV( LDT, LDT ), QTMP( LDT, LDT ), T( LDT, LDT ), TMP( LDT, LDT ), TSAV( LDT, LDT ), TSAV1( LDT, LDT ), TTMP( LDT, LDT ), W( LDT ), WORK( LWORK ), WTMP( LDT );
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEHRD, CHSEQR, CHST01, CLACPY, CSSCAL, CTRSEN, CUNGHR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, MAX, REAL, SQRT

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // EPSIN = 2**(-24) = precision to which input data computed

      EPS = max( EPS, EPSIN );
      RMAX[1] = ZERO;
      RMAX[2] = ZERO;
      RMAX[3] = ZERO;
      LMAX[1] = 0;
      LMAX[2] = 0;
      LMAX[3] = 0;
      KNT = 0;
      NINFO[1] = 0;
      NINFO[2] = 0;
      NINFO[3] = 0;
      VAL[1] = sqrt( SMLNUM );
      VAL[2] = ONE;
      VAL[3] = sqrt( sqrt( BIGNUM ) );

      // Read input data until N=0.  Assume input eigenvalues are sorted
      // lexicographically (increasing by real part, then decreasing by
      // imaginary part)

      } // 10
      READ( NIN, FMT = * )N, NDIM, ISRT;
      if (N == 0) return;
      READ( NIN, FMT = * )( ISELEC( I ), I = 1, NDIM );
      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N );
      } // 20
      READ( NIN, FMT = * )SIN, SEPIN;

      TNRM = CLANGE( 'M', N, N, TMP, LDT, RWORK );
      for (ISCL = 1; ISCL <= 3; ISCL++) { // 200

         // Scale input matrix

         KNT = KNT + 1;
         clacpy('F', N, N, TMP, LDT, T, LDT );
         VMUL = VAL( ISCL );
         for (I = 1; I <= N; I++) { // 30
            csscal(N, VMUL, T( 1, I ), 1 );
         } // 30
         if (TNRM == ZERO) VMUL = ONE;
         clacpy('F', N, N, T, LDT, TSAV, LDT );

         // Compute Schur form

         cgehrd(N, 1, N, T, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, INFO );
         if ( INFO != 0 ) {
            LMAX[1] = KNT;
            NINFO[1] = NINFO( 1 ) + 1;
            GO TO 200;
         }

         // Generate unitary matrix

         clacpy('L', N, N, T, LDT, Q, LDT );
         cunghr(N, 1, N, Q, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, INFO );

         // Compute Schur form

         for (J = 1; J <= N - 2; J++) { // 50
            for (I = J + 2; I <= N; I++) { // 40
               T[I][J] = CZERO;
            } // 40
         } // 50
         chseqr('S', 'V', N, 1, N, T, LDT, W, Q, LDT, WORK, LWORK, INFO );
         if ( INFO != 0 ) {
            LMAX[2] = KNT;
            NINFO[2] = NINFO( 2 ) + 1;
            GO TO 200;
         }

         // Sort, select eigenvalues

         for (I = 1; I <= N; I++) { // 60
            IPNT[I] = I;
            SELECT[I] = false;
         } // 60
         if ( ISRT == 0 ) {
            for (I = 1; I <= N; I++) { // 70
               WSRT[I] = double( W( I ) );
            } // 70
         } else {
            for (I = 1; I <= N; I++) { // 80
               WSRT[I] = AIMAG( W( I ) );
            } // 80
         }
         for (I = 1; I <= N - 1; I++) { // 100
            KMIN = I;
            VMIN = WSRT( I );
            for (J = I + 1; J <= N; J++) { // 90
               if ( WSRT( J ) < VMIN ) {
                  KMIN = J;
                  VMIN = WSRT( J );
               }
            } // 90
            WSRT[KMIN] = WSRT( I );
            WSRT[I] = VMIN;
            ITMP = IPNT( I );
            IPNT[I] = IPNT( KMIN );
            IPNT[KMIN] = ITMP;
         } // 100
         for (I = 1; I <= NDIM; I++) { // 110
            SELECT[IPNT( ISELEC( I ) )] = true;
         } // 110

         // Compute condition numbers

         clacpy('F', N, N, Q, LDT, QSAV, LDT );
         clacpy('F', N, N, T, LDT, TSAV1, LDT );
         ctrsen('B', 'V', SELECT, N, T, LDT, Q, LDT, WTMP, M, S, SEP, WORK, LWORK, INFO );
         if ( INFO != 0 ) {
            LMAX[3] = KNT;
            NINFO[3] = NINFO( 3 ) + 1;
            GO TO 200;
         }
         SEPTMP = SEP / VMUL;
         STMP = S;

         // Compute residuals

         chst01(N, 1, N, TSAV, LDT, T, LDT, Q, LDT, WORK, LWORK, RWORK, RESULT );
         VMAX = max( RESULT( 1 ), RESULT( 2 ) );
         if ( VMAX > RMAX( 1 ) ) {
            RMAX[1] = VMAX;
            if( NINFO( 1 ) == 0 ) LMAX( 1 ) = KNT;
         }

         // Compare condition number for eigenvalue cluster
         // taking its condition number into account

         V = max( TWO*double( N )*EPS*TNRM, SMLNUM );
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
         TOL = max( TOL, SMLNUM / EPS );
         TOLIN = max( TOLIN, SMLNUM / EPS );
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
            RMAX[2] = VMAX;
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
         TOL = max( TOL, SMLNUM / EPS );
         TOLIN = max( TOLIN, SMLNUM / EPS );
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
            RMAX[2] = VMAX;
            if( NINFO( 2 ) == 0 ) LMAX( 2 ) = KNT;
         }

         // Compare condition number for eigenvalue cluster
         // without taking its condition number into account

         if ( SIN <= REAL( 2*N )*EPS && STMP <= double( 2*N )*EPS ) {
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
            RMAX[3] = VMAX;
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
            RMAX[3] = VMAX;
            if( NINFO( 3 ) == 0 ) LMAX( 3 ) = KNT;
         }

         // Compute eigenvalue condition number only and compare
         // Update Q

         VMAX = ZERO;
         clacpy('F', N, N, TSAV1, LDT, TTMP, LDT );
         clacpy('F', N, N, QSAV, LDT, QTMP, LDT );
         SEPTMP = -ONE;
         STMP = -ONE;
         ctrsen('E', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WTMP, M, STMP, SEPTMP, WORK, LWORK, INFO );
         if ( INFO != 0 ) {
            LMAX[3] = KNT;
            NINFO[3] = NINFO( 3 ) + 1;
            GO TO 200;
         }
         if (S != STMP) VMAX = ONE / EPS;
         IF( -ONE != SEPTMP ) VMAX = ONE / EPS;
         for (I = 1; I <= N; I++) { // 130
            for (J = 1; J <= N; J++) { // 120
               if( TTMP( I, J ) != T( I, J ) ) VMAX = ONE / EPS;
               IF( QTMP( I, J ) != Q( I, J ) ) VMAX = ONE / EPS;
            } // 120
         } // 130

         // Compute invariant subspace condition number only and compare
         // Update Q

         clacpy('F', N, N, TSAV1, LDT, TTMP, LDT );
         clacpy('F', N, N, QSAV, LDT, QTMP, LDT );
         SEPTMP = -ONE;
         STMP = -ONE;
         ctrsen('V', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WTMP, M, STMP, SEPTMP, WORK, LWORK, INFO );
         if ( INFO != 0 ) {
            LMAX[3] = KNT;
            NINFO[3] = NINFO( 3 ) + 1;
            GO TO 200;
         }
         if (-ONE != STMP) VMAX = ONE / EPS;
         IF( SEP != SEPTMP ) VMAX = ONE / EPS;
         for (I = 1; I <= N; I++) { // 150
            for (J = 1; J <= N; J++) { // 140
               if( TTMP( I, J ) != T( I, J ) ) VMAX = ONE / EPS;
               IF( QTMP( I, J ) != Q( I, J ) ) VMAX = ONE / EPS;
            } // 140
         } // 150

         // Compute eigenvalue condition number only and compare
         // Do not update Q

         clacpy('F', N, N, TSAV1, LDT, TTMP, LDT );
         clacpy('F', N, N, QSAV, LDT, QTMP, LDT );
         SEPTMP = -ONE;
         STMP = -ONE;
         ctrsen('E', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WTMP, M, STMP, SEPTMP, WORK, LWORK, INFO );
         if ( INFO != 0 ) {
            LMAX[3] = KNT;
            NINFO[3] = NINFO( 3 ) + 1;
            GO TO 200;
         }
         if (S != STMP) VMAX = ONE / EPS;
         IF( -ONE != SEPTMP ) VMAX = ONE / EPS;
         for (I = 1; I <= N; I++) { // 170
            for (J = 1; J <= N; J++) { // 160
               if( TTMP( I, J ) != T( I, J ) ) VMAX = ONE / EPS;
               IF( QTMP( I, J ) != QSAV( I, J ) ) VMAX = ONE / EPS;
            } // 160
         } // 170

         // Compute invariant subspace condition number only and compare
         // Do not update Q

         clacpy('F', N, N, TSAV1, LDT, TTMP, LDT );
         clacpy('F', N, N, QSAV, LDT, QTMP, LDT );
         SEPTMP = -ONE;
         STMP = -ONE;
         ctrsen('V', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WTMP, M, STMP, SEPTMP, WORK, LWORK, INFO );
         if ( INFO != 0 ) {
            LMAX[3] = KNT;
            NINFO[3] = NINFO( 3 ) + 1;
            GO TO 200;
         }
         if (-ONE != STMP) VMAX = ONE / EPS;
         IF( SEP != SEPTMP ) VMAX = ONE / EPS;
         for (I = 1; I <= N; I++) { // 190
            for (J = 1; J <= N; J++) { // 180
               if( TTMP( I, J ) != T( I, J ) ) VMAX = ONE / EPS;
               IF( QTMP( I, J ) != QSAV( I, J ) ) VMAX = ONE / EPS;
            } // 180
         } // 190
         if ( VMAX > RMAX( 1 ) ) {
            RMAX[1] = VMAX;
            if( NINFO( 1 ) == 0 ) LMAX( 1 ) = KNT;
         }
      } // 200
      GO TO 10;
      }
