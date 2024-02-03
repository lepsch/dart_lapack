      SUBROUTINE ZGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOB;
      int                IHI, ILO, INFO, LDA, LDB, N;
      // ..
      // .. Array Arguments ..
      double             LSCALE( * ), RSCALE( * ), WORK( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, HALF, ONE;
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0 ;
      double             THREE, SCLFAC;
      const              THREE = 3.0, SCLFAC = 1.0e+1 ;
      COMPLEX*16         CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, ICAB, IFLOW, IP1, IR, IRAB, IT, J, JC, JP1, K, KOUNT, L, LCAB, LM1, LRAB, LSFMAX, LSFMIN, M, NR, NRP2;
      double             ALPHA, BASL, BETA, CAB, CMAX, COEF, COEF2, COEF5, COR, EW, EWC, GAMMA, PGAMMA, RAB, SFMAX, SFMIN, SUM, T, TA, TB, TC;
      COMPLEX*16         CDUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DDOT, DLAMCH;
      // EXTERNAL LSAME, IZAMAX, DDOT, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DSCAL, XERBLA, ZDSCAL, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, INT, LOG10, MAX, MIN, SIGN
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) );
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0;
      if ( !LSAME( JOB, 'N' ) && !LSAME( JOB, 'P' ) && !LSAME( JOB, 'S' ) && !LSAME( JOB, 'B' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4;
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('ZGGBAL', -INFO );
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         ILO = 1;
         IHI = N;
         return;
      }

      if ( N == 1 ) {
         ILO = 1;
         IHI = N;
         LSCALE( 1 ) = ONE;
         RSCALE( 1 ) = ONE;
         return;
      }

      if ( LSAME( JOB, 'N' ) ) {
         ILO = 1;
         IHI = N;
         for (I = 1; I <= N; I++) { // 10
            LSCALE( I ) = ONE;
            RSCALE( I ) = ONE;
         } // 10
         return;
      }

      K = 1;
      L = N;
      IF( LSAME( JOB, 'S' ) ) GO TO 190;

      GO TO 30;

      // Permute the matrices A and B to isolate the eigenvalues.

      // Find row with one nonzero in columns 1 through L

      } // 20
      L = LM1;
      if (L != 1) GO TO 30;

      RSCALE( 1 ) = 1;
      LSCALE( 1 ) = 1;
      GO TO 190;

      } // 30
      LM1 = L - 1;
      DO 80 I = L, 1, -1;
         for (J = 1; J <= LM1; J++) { // 40
            JP1 = J + 1;
            IF( A( I, J ) != CZERO || B( I, J ) != CZERO ) GO TO 50;
         } // 40
         J = L;
         GO TO 70;

         } // 50
         for (J = JP1; J <= L; J++) { // 60
            IF( A( I, J ) != CZERO || B( I, J ) != CZERO ) GO TO 80;
         } // 60
         J = JP1 - 1;

         } // 70
         M = L;
         IFLOW = 1;
         GO TO 160;
      } // 80
      GO TO 100;

      // Find column with one nonzero in rows K through N

      } // 90
      K = K + 1;

      } // 100
      for (J = K; J <= L; J++) { // 150
         for (I = K; I <= LM1; I++) { // 110
            IP1 = I + 1;
            IF( A( I, J ) != CZERO || B( I, J ) != CZERO ) GO TO 120;
         } // 110
         I = L;
         GO TO 140;
         } // 120
         for (I = IP1; I <= L; I++) { // 130
            IF( A( I, J ) != CZERO || B( I, J ) != CZERO ) GO TO 150;
         } // 130
         I = IP1 - 1;
         } // 140
         M = K;
         IFLOW = 2;
         GO TO 160;
      } // 150
      GO TO 190;

      // Permute rows M and I

      } // 160
      LSCALE( M ) = I;
      if (I == M) GO TO 170;
      zswap(N-K+1, A( I, K ), LDA, A( M, K ), LDA );
      zswap(N-K+1, B( I, K ), LDB, B( M, K ), LDB );

      // Permute columns M and J

      } // 170
      RSCALE( M ) = J;
      if (J == M) GO TO 180;
      zswap(L, A( 1, J ), 1, A( 1, M ), 1 );
      zswap(L, B( 1, J ), 1, B( 1, M ), 1 );

      } // 180
      GO TO ( 20, 90 )IFLOW;

      } // 190
      ILO = K;
      IHI = L;

      if ( LSAME( JOB, 'P' ) ) {
         for (I = ILO; I <= IHI; I++) { // 195
            LSCALE( I ) = ONE;
            RSCALE( I ) = ONE;
         } // 195
         return;
      }

      if (ILO == IHI) RETURN;

      // Balance the submatrix in rows ILO to IHI.

      NR = IHI - ILO + 1;
      for (I = ILO; I <= IHI; I++) { // 200
         RSCALE( I ) = ZERO;
         LSCALE( I ) = ZERO;

         WORK( I ) = ZERO;
         WORK( I+N ) = ZERO;
         WORK( I+2*N ) = ZERO;
         WORK( I+3*N ) = ZERO;
         WORK( I+4*N ) = ZERO;
         WORK( I+5*N ) = ZERO;
      } // 200

      // Compute right side vector in resulting linear equations

      BASL = LOG10( SCLFAC );
      for (I = ILO; I <= IHI; I++) { // 240
         for (J = ILO; J <= IHI; J++) { // 230
            if ( A( I, J ) == CZERO ) {
               TA = ZERO;
               GO TO 210;
            }
            TA = LOG10( CABS1( A( I, J ) ) ) / BASL;

            } // 210
            if ( B( I, J ) == CZERO ) {
               TB = ZERO;
               GO TO 220;
            }
            TB = LOG10( CABS1( B( I, J ) ) ) / BASL;

            } // 220
            WORK( I+4*N ) = WORK( I+4*N ) - TA - TB;
            WORK( J+5*N ) = WORK( J+5*N ) - TA - TB;
         } // 230
      } // 240

      COEF = ONE / DBLE( 2*NR );
      COEF2 = COEF*COEF;
      COEF5 = HALF*COEF2;
      NRP2 = NR + 2;
      BETA = ZERO;
      IT = 1;

      // Start generalized conjugate gradient iteration

      } // 250

      GAMMA = DDOT( NR, WORK( ILO+4*N ), 1, WORK( ILO+4*N ), 1 ) + DDOT( NR, WORK( ILO+5*N ), 1, WORK( ILO+5*N ), 1 );

      EW = ZERO;
      EWC = ZERO;
      for (I = ILO; I <= IHI; I++) { // 260
         EW = EW + WORK( I+4*N );
         EWC = EWC + WORK( I+5*N );
      } // 260

      GAMMA = COEF*GAMMA - COEF2*( EW**2+EWC**2 ) - COEF5*( EW-EWC )**2;
      if (GAMMA == ZERO) GO TO 350       IF( IT != 1 ) BETA = GAMMA / PGAMMA;
      T = COEF5*( EWC-THREE*EW );
      TC = COEF5*( EW-THREE*EWC );

      dscal(NR, BETA, WORK( ILO ), 1 );
      dscal(NR, BETA, WORK( ILO+N ), 1 );

      daxpy(NR, COEF, WORK( ILO+4*N ), 1, WORK( ILO+N ), 1 );
      daxpy(NR, COEF, WORK( ILO+5*N ), 1, WORK( ILO ), 1 );

      for (I = ILO; I <= IHI; I++) { // 270
         WORK( I ) = WORK( I ) + TC;
         WORK( I+N ) = WORK( I+N ) + T;
      } // 270

      // Apply matrix to vector

      for (I = ILO; I <= IHI; I++) { // 300
         KOUNT = 0;
         SUM = ZERO;
         for (J = ILO; J <= IHI; J++) { // 290
            IF( A( I, J ) == CZERO ) GO TO 280;
            KOUNT = KOUNT + 1;
            SUM = SUM + WORK( J );
            } // 280
            IF( B( I, J ) == CZERO ) GO TO 290;
            KOUNT = KOUNT + 1;
            SUM = SUM + WORK( J );
         } // 290
         WORK( I+2*N ) = DBLE( KOUNT )*WORK( I+N ) + SUM;
      } // 300

      for (J = ILO; J <= IHI; J++) { // 330
         KOUNT = 0;
         SUM = ZERO;
         for (I = ILO; I <= IHI; I++) { // 320
            IF( A( I, J ) == CZERO ) GO TO 310;
            KOUNT = KOUNT + 1;
            SUM = SUM + WORK( I+N );
            } // 310
            IF( B( I, J ) == CZERO ) GO TO 320;
            KOUNT = KOUNT + 1;
            SUM = SUM + WORK( I+N );
         } // 320
         WORK( J+3*N ) = DBLE( KOUNT )*WORK( J ) + SUM;
      } // 330

      SUM = DDOT( NR, WORK( ILO+N ), 1, WORK( ILO+2*N ), 1 ) + DDOT( NR, WORK( ILO ), 1, WORK( ILO+3*N ), 1 );
      ALPHA = GAMMA / SUM;

      // Determine correction to current iteration

      CMAX = ZERO;
      for (I = ILO; I <= IHI; I++) { // 340
         COR = ALPHA*WORK( I+N );
         IF( ABS( COR ) > CMAX ) CMAX = ABS( COR );
         LSCALE( I ) = LSCALE( I ) + COR;
         COR = ALPHA*WORK( I );
         IF( ABS( COR ) > CMAX ) CMAX = ABS( COR );
         RSCALE( I ) = RSCALE( I ) + COR;
      } // 340
      if (CMAX < HALF) GO TO 350;

      daxpy(NR, -ALPHA, WORK( ILO+2*N ), 1, WORK( ILO+4*N ), 1 );
      daxpy(NR, -ALPHA, WORK( ILO+3*N ), 1, WORK( ILO+5*N ), 1 );

      PGAMMA = GAMMA;
      IT = IT + 1;
      if (IT <= NRP2) GO TO 250;

      // End generalized conjugate gradient iteration

      } // 350
      SFMIN = DLAMCH( 'S' );
      SFMAX = ONE / SFMIN;
      LSFMIN = INT( LOG10( SFMIN ) / BASL+ONE );
      LSFMAX = INT( LOG10( SFMAX ) / BASL );
      for (I = ILO; I <= IHI; I++) { // 360
         IRAB = IZAMAX( N-ILO+1, A( I, ILO ), LDA );
         RAB = ABS( A( I, IRAB+ILO-1 ) );
         IRAB = IZAMAX( N-ILO+1, B( I, ILO ), LDB );
         RAB = MAX( RAB, ABS( B( I, IRAB+ILO-1 ) ) );
         LRAB = INT( LOG10( RAB+SFMIN ) / BASL+ONE );
         IR = INT(LSCALE( I ) + SIGN( HALF, LSCALE( I ) ));
         IR = MIN( MAX( IR, LSFMIN ), LSFMAX, LSFMAX-LRAB );
         LSCALE( I ) = SCLFAC**IR;
         ICAB = IZAMAX( IHI, A( 1, I ), 1 );
         CAB = ABS( A( ICAB, I ) );
         ICAB = IZAMAX( IHI, B( 1, I ), 1 );
         CAB = MAX( CAB, ABS( B( ICAB, I ) ) );
         LCAB = INT( LOG10( CAB+SFMIN ) / BASL+ONE );
         JC = INT(RSCALE( I ) + SIGN( HALF, RSCALE( I ) ));
         JC = MIN( MAX( JC, LSFMIN ), LSFMAX, LSFMAX-LCAB );
         RSCALE( I ) = SCLFAC**JC;
      } // 360

      // Row scaling of matrices A and B

      for (I = ILO; I <= IHI; I++) { // 370
         zdscal(N-ILO+1, LSCALE( I ), A( I, ILO ), LDA );
         zdscal(N-ILO+1, LSCALE( I ), B( I, ILO ), LDB );
      } // 370

      // Column scaling of matrices A and B

      for (J = ILO; J <= IHI; J++) { // 380
         zdscal(IHI, RSCALE( J ), A( 1, J ), 1 );
         zdscal(IHI, RSCALE( J ), B( 1, J ), 1 );
      } // 380

      return;

      // End of ZGGBAL

      }
