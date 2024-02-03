      SUBROUTINE SGET38( RMAX, LMAX, NINFO, KNT, NIN )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, NIN;
      // ..
      // .. Array Arguments ..
      int                LMAX( 3 ), NINFO( 3 );
      REAL               RMAX( 3 )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO
      const              ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 ;
      REAL               EPSIN
      const              EPSIN = 5.9605E-8 ;
      int                LDT, LWORK;
      const              LDT = 20, LWORK = 2*LDT*( 10+LDT ) ;
      int                LIWORK;
      const              LIWORK = LDT*LDT ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, ISCL, ITMP, J, KMIN, M, N, NDIM;
      REAL               BIGNUM, EPS, S, SEP, SEPIN, SEPTMP, SIN, SMLNUM, STMP, TNRM, TOL, TOLIN, V, VIMIN, VMAX, VMUL, VRMIN
      // ..
      // .. Local Arrays ..
      bool               SELECT( LDT );
      int                IPNT( LDT ), ISELEC( LDT ), IWORK( LIWORK );
      REAL               Q( LDT, LDT ), QSAV( LDT, LDT ), QTMP( LDT, LDT ), RESULT( 2 ), T( LDT, LDT ), TMP( LDT, LDT ), TSAV( LDT, LDT ), TSAV1( LDT, LDT ), TTMP( LDT, LDT ), VAL( 3 ), WI( LDT ), WITMP( LDT ), WORK( LWORK ), WR( LDT ), WRTMP( LDT )
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEHRD, SHSEQR, SHST01, SLACPY, SORGHR, SSCAL, STRSEN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL, SQRT
      // ..
      // .. Executable Statements ..

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM

      // EPSIN = 2**(-24) = precision to which input data computed

      EPS = MAX( EPS, EPSIN )
      RMAX( 1 ) = ZERO
      RMAX( 2 ) = ZERO
      RMAX( 3 ) = ZERO
      LMAX( 1 ) = 0
      LMAX( 2 ) = 0
      LMAX( 3 ) = 0
      KNT = 0
      NINFO( 1 ) = 0
      NINFO( 2 ) = 0
      NINFO( 3 ) = 0

      VAL( 1 ) = SQRT( SMLNUM )
      VAL( 2 ) = ONE
      VAL( 3 ) = SQRT( SQRT( BIGNUM ) )

      // Read input data until N=0.  Assume input eigenvalues are sorted
      // lexicographically (increasing by real part, then decreasing by
      // imaginary part)

   10 CONTINUE
      READ( NIN, FMT = * )N, NDIM
      IF( N.EQ.0 ) RETURN
      READ( NIN, FMT = * )( ISELEC( I ), I = 1, NDIM )
      DO 20 I = 1, N
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N )
   20 CONTINUE
      READ( NIN, FMT = * )SIN, SEPIN

      TNRM = SLANGE( 'M', N, N, TMP, LDT, WORK )
      DO 160 ISCL = 1, 3

         // Scale input matrix

         KNT = KNT + 1
         CALL SLACPY( 'F', N, N, TMP, LDT, T, LDT )
         VMUL = VAL( ISCL )
         DO 30 I = 1, N
            CALL SSCAL( N, VMUL, T( 1, I ), 1 )
   30    CONTINUE
         IF( TNRM.EQ.ZERO ) VMUL = ONE
         CALL SLACPY( 'F', N, N, T, LDT, TSAV, LDT )

         // Compute Schur form

         CALL SGEHRD( N, 1, N, T, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 1 ) = KNT
            NINFO( 1 ) = NINFO( 1 ) + 1
            GO TO 160
         }

         // Generate orthogonal matrix

         CALL SLACPY( 'L', N, N, T, LDT, Q, LDT )
         CALL SORGHR( N, 1, N, Q, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, INFO )

         // Compute Schur form

         CALL SHSEQR( 'S', 'V', N, 1, N, T, LDT, WR, WI, Q, LDT, WORK, LWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 2 ) = KNT
            NINFO( 2 ) = NINFO( 2 ) + 1
            GO TO 160
         }

         // Sort, select eigenvalues

         DO 40 I = 1, N
            IPNT( I ) = I
            SELECT( I ) = .FALSE.
   40    CONTINUE
         CALL SCOPY( N, WR, 1, WRTMP, 1 )
         CALL SCOPY( N, WI, 1, WITMP, 1 )
         DO 60 I = 1, N - 1
            KMIN = I
            VRMIN = WRTMP( I )
            VIMIN = WITMP( I )
            DO 50 J = I + 1, N
               if ( WRTMP( J ).LT.VRMIN ) {
                  KMIN = J
                  VRMIN = WRTMP( J )
                  VIMIN = WITMP( J )
               }
   50       CONTINUE
            WRTMP( KMIN ) = WRTMP( I )
            WITMP( KMIN ) = WITMP( I )
            WRTMP( I ) = VRMIN
            WITMP( I ) = VIMIN
            ITMP = IPNT( I )
            IPNT( I ) = IPNT( KMIN )
            IPNT( KMIN ) = ITMP
   60    CONTINUE
         DO 70 I = 1, NDIM
            SELECT( IPNT( ISELEC( I ) ) ) = .TRUE.
   70    CONTINUE

         // Compute condition numbers

         CALL SLACPY( 'F', N, N, Q, LDT, QSAV, LDT )
         CALL SLACPY( 'F', N, N, T, LDT, TSAV1, LDT )
         CALL STRSEN( 'B', 'V', SELECT, N, T, LDT, Q, LDT, WRTMP, WITMP, M, S, SEP, WORK, LWORK, IWORK, LIWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 160
         }
         SEPTMP = SEP / VMUL
         STMP = S

         // Compute residuals

         CALL SHST01( N, 1, N, TSAV, LDT, T, LDT, Q, LDT, WORK, LWORK, RESULT )
         VMAX = MAX( RESULT( 1 ), RESULT( 2 ) )
         if ( VMAX.GT.RMAX( 1 ) ) {
            RMAX( 1 ) = VMAX
            IF( NINFO( 1 ).EQ.0 ) LMAX( 1 ) = KNT
         }

         // Compare condition number for eigenvalue cluster
        t // aking its condition number into account

         V = MAX( TWO*REAL( N )*EPS*TNRM, SMLNUM )
         IF( TNRM.EQ.ZERO ) V = ONE
         if ( V.GT.SEPTMP ) {
            TOL = ONE
         } else {
            TOL = V / SEPTMP
         }
         if ( V.GT.SEPIN ) {
            TOLIN = ONE
         } else {
            TOLIN = V / SEPIN
         }
         TOL = MAX( TOL, SMLNUM / EPS )
         TOLIN = MAX( TOLIN, SMLNUM / EPS )
         if ( EPS*( SIN-TOLIN ).GT.STMP+TOL ) {
            VMAX = ONE / EPS
         } else if ( SIN-TOLIN.GT.STMP+TOL ) {
            VMAX = ( SIN-TOLIN ) / ( STMP+TOL )
         } else if ( SIN+TOLIN.LT.EPS*( STMP-TOL ) ) {
            VMAX = ONE / EPS
         } else if ( SIN+TOLIN.LT.STMP-TOL ) {
            VMAX = ( STMP-TOL ) / ( SIN+TOLIN )
         } else {
            VMAX = ONE
         }
         if ( VMAX.GT.RMAX( 2 ) ) {
            RMAX( 2 ) = VMAX
            IF( NINFO( 2 ).EQ.0 ) LMAX( 2 ) = KNT
         }

         // Compare condition numbers for invariant subspace
        t // aking its condition number into account

         if ( V.GT.SEPTMP*STMP ) {
            TOL = SEPTMP
         } else {
            TOL = V / STMP
         }
         if ( V.GT.SEPIN*SIN ) {
            TOLIN = SEPIN
         } else {
            TOLIN = V / SIN
         }
         TOL = MAX( TOL, SMLNUM / EPS )
         TOLIN = MAX( TOLIN, SMLNUM / EPS )
         if ( EPS*( SEPIN-TOLIN ).GT.SEPTMP+TOL ) {
            VMAX = ONE / EPS
         } else if ( SEPIN-TOLIN.GT.SEPTMP+TOL ) {
            VMAX = ( SEPIN-TOLIN ) / ( SEPTMP+TOL )
         } else if ( SEPIN+TOLIN.LT.EPS*( SEPTMP-TOL ) ) {
            VMAX = ONE / EPS
         } else if ( SEPIN+TOLIN.LT.SEPTMP-TOL ) {
            VMAX = ( SEPTMP-TOL ) / ( SEPIN+TOLIN )
         } else {
            VMAX = ONE
         }
         if ( VMAX.GT.RMAX( 2 ) ) {
            RMAX( 2 ) = VMAX
            IF( NINFO( 2 ).EQ.0 ) LMAX( 2 ) = KNT
         }

         // Compare condition number for eigenvalue cluster
         // without taking its condition number into account

         if ( SIN.LE.REAL( 2*N )*EPS .AND. STMP.LE.REAL( 2*N )*EPS ) {
            VMAX = ONE
         } else if ( EPS*SIN.GT.STMP ) {
            VMAX = ONE / EPS
         } else if ( SIN.GT.STMP ) {
            VMAX = SIN / STMP
         } else if ( SIN.LT.EPS*STMP ) {
            VMAX = ONE / EPS
         } else if ( SIN.LT.STMP ) {
            VMAX = STMP / SIN
         } else {
            VMAX = ONE
         }
         if ( VMAX.GT.RMAX( 3 ) ) {
            RMAX( 3 ) = VMAX
            IF( NINFO( 3 ).EQ.0 ) LMAX( 3 ) = KNT
         }

         // Compare condition numbers for invariant subspace
         // without taking its condition number into account

         if ( SEPIN.LE.V .AND. SEPTMP.LE.V ) {
            VMAX = ONE
         } else if ( EPS*SEPIN.GT.SEPTMP ) {
            VMAX = ONE / EPS
         } else if ( SEPIN.GT.SEPTMP ) {
            VMAX = SEPIN / SEPTMP
         } else if ( SEPIN.LT.EPS*SEPTMP ) {
            VMAX = ONE / EPS
         } else if ( SEPIN.LT.SEPTMP ) {
            VMAX = SEPTMP / SEPIN
         } else {
            VMAX = ONE
         }
         if ( VMAX.GT.RMAX( 3 ) ) {
            RMAX( 3 ) = VMAX
            IF( NINFO( 3 ).EQ.0 ) LMAX( 3 ) = KNT
         }

         // Compute eigenvalue condition number only and compare
         // Update Q

         VMAX = ZERO
         CALL SLACPY( 'F', N, N, TSAV1, LDT, TTMP, LDT )
         CALL SLACPY( 'F', N, N, QSAV, LDT, QTMP, LDT )
         SEPTMP = -ONE
         STMP = -ONE
         CALL STRSEN( 'E', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, LIWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 160
         }
         IF( S.NE.STMP ) VMAX = ONE / EPS          IF( -ONE.NE.SEPTMP ) VMAX = ONE / EPS
         DO 90 I = 1, N
            DO 80 J = 1, N
               IF( TTMP( I, J ).NE.T( I, J ) ) VMAX = ONE / EPS                IF( QTMP( I, J ).NE.Q( I, J ) ) VMAX = ONE / EPS
   80       CONTINUE
   90    CONTINUE

         // Compute invariant subspace condition number only and compare
         // Update Q

         CALL SLACPY( 'F', N, N, TSAV1, LDT, TTMP, LDT )
         CALL SLACPY( 'F', N, N, QSAV, LDT, QTMP, LDT )
         SEPTMP = -ONE
         STMP = -ONE
         CALL STRSEN( 'V', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, LIWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 160
         }
         IF( -ONE.NE.STMP ) VMAX = ONE / EPS          IF( SEP.NE.SEPTMP ) VMAX = ONE / EPS
         DO 110 I = 1, N
            DO 100 J = 1, N
               IF( TTMP( I, J ).NE.T( I, J ) ) VMAX = ONE / EPS                IF( QTMP( I, J ).NE.Q( I, J ) ) VMAX = ONE / EPS
  100       CONTINUE
  110    CONTINUE

         // Compute eigenvalue condition number only and compare
         // Do not update Q

         CALL SLACPY( 'F', N, N, TSAV1, LDT, TTMP, LDT )
         CALL SLACPY( 'F', N, N, QSAV, LDT, QTMP, LDT )
         SEPTMP = -ONE
         STMP = -ONE
         CALL STRSEN( 'E', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, LIWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 160
         }
         IF( S.NE.STMP ) VMAX = ONE / EPS          IF( -ONE.NE.SEPTMP ) VMAX = ONE / EPS
         DO 130 I = 1, N
            DO 120 J = 1, N
               IF( TTMP( I, J ).NE.T( I, J ) ) VMAX = ONE / EPS                IF( QTMP( I, J ).NE.QSAV( I, J ) ) VMAX = ONE / EPS
  120       CONTINUE
  130    CONTINUE

         // Compute invariant subspace condition number only and compare
         // Do not update Q

         CALL SLACPY( 'F', N, N, TSAV1, LDT, TTMP, LDT )
         CALL SLACPY( 'F', N, N, QSAV, LDT, QTMP, LDT )
         SEPTMP = -ONE
         STMP = -ONE
         CALL STRSEN( 'V', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WRTMP, WITMP, M, STMP, SEPTMP, WORK, LWORK, IWORK, LIWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 160
         }
         IF( -ONE.NE.STMP ) VMAX = ONE / EPS          IF( SEP.NE.SEPTMP ) VMAX = ONE / EPS
         DO 150 I = 1, N
            DO 140 J = 1, N
               IF( TTMP( I, J ).NE.T( I, J ) ) VMAX = ONE / EPS                IF( QTMP( I, J ).NE.QSAV( I, J ) ) VMAX = ONE / EPS
  140       CONTINUE
  150    CONTINUE
         if ( VMAX.GT.RMAX( 1 ) ) {
            RMAX( 1 ) = VMAX
            IF( NINFO( 1 ).EQ.0 ) LMAX( 1 ) = KNT
         }
  160 CONTINUE
      GO TO 10

      // End of SGET38

      }
