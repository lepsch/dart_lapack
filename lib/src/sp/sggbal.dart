      void sggbal(JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOB;
      int                IHI, ILO, INFO, LDA, LDB, N;
      double               A( LDA, * ), B( LDB, * ), LSCALE( * ), RSCALE( * ), WORK( * );
      // ..

      double               ZERO, HALF, ONE;
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0 ;
      double               THREE, SCLFAC;
      const              THREE = 3.0, SCLFAC = 1.0e+1 ;
      int                I, ICAB, IFLOW, IP1, IR, IRAB, IT, J, JC, JP1, K, KOUNT, L, LCAB, LM1, LRAB, LSFMAX, LSFMIN, M, NR, NRP2;
      double               ALPHA, BASL, BETA, CAB, CMAX, COEF, COEF2, COEF5, COR, EW, EWC, GAMMA, PGAMMA, RAB, SFMAX, SFMIN, SUM, T, TA, TB, TC;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ISAMAX;
      //- REAL               SDOT, SLAMCH;
      // EXTERNAL lsame, ISAMAX, SDOT, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SSCAL, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, LOG10, MAX, MIN, REAL, SIGN

      // Test the input parameters

      INFO = 0;
      if ( !lsame( JOB, 'N' ) && !lsame( JOB, 'P' ) && !lsame( JOB, 'S' ) && !lsame( JOB, 'B' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('SGGBAL', -INFO );
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
         LSCALE[1] = ONE;
         RSCALE[1] = ONE;
         return;
      }

      if ( lsame( JOB, 'N' ) ) {
         ILO = 1;
         IHI = N;
         for (I = 1; I <= N; I++) { // 10
            LSCALE[I] = ONE;
            RSCALE[I] = ONE;
         } // 10
         return;
      }

      K = 1;
      L = N;
      if( lsame( JOB, 'S' ) ) GO TO 190;

      GO TO 30;

      // Permute the matrices A and B to isolate the eigenvalues.

      // Find row with one nonzero in columns 1 through L

      } // 20
      L = LM1;
      if (L != 1) GO TO 30;

      RSCALE[1] = ONE;
      LSCALE[1] = ONE;
      GO TO 190;

      } // 30
      LM1 = L - 1;
      for (I = L; I >= 1; I--) { // 80
         for (J = 1; J <= LM1; J++) { // 40
            JP1 = J + 1;
            if( A( I, J ) != ZERO || B( I, J ) != ZERO ) GO TO 50;
         } // 40
         J = L;
         GO TO 70;

         } // 50
         for (J = JP1; J <= L; J++) { // 60
            if( A( I, J ) != ZERO || B( I, J ) != ZERO ) GO TO 80;
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
            if( A( I, J ) != ZERO || B( I, J ) != ZERO ) GO TO 120;
         } // 110
         I = L;
         GO TO 140;
         } // 120
         for (I = IP1; I <= L; I++) { // 130
            if( A( I, J ) != ZERO || B( I, J ) != ZERO ) GO TO 150;
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
      LSCALE[M] = I;
      if (I == M) GO TO 170;
      sswap(N-K+1, A( I, K ), LDA, A( M, K ), LDA );
      sswap(N-K+1, B( I, K ), LDB, B( M, K ), LDB );

      // Permute columns M and J

      } // 170
      RSCALE[M] = J;
      if (J == M) GO TO 180;
      sswap(L, A( 1, J ), 1, A( 1, M ), 1 );
      sswap(L, B( 1, J ), 1, B( 1, M ), 1 );

      } // 180
      GO TO ( 20, 90 )IFLOW;

      } // 190
      ILO = K;
      IHI = L;

      if ( lsame( JOB, 'P' ) ) {
         for (I = ILO; I <= IHI; I++) { // 195
            LSCALE[I] = ONE;
            RSCALE[I] = ONE;
         } // 195
         return;
      }

      if (ILO == IHI) return;

      // Balance the submatrix in rows ILO to IHI.

      NR = IHI - ILO + 1;
      for (I = ILO; I <= IHI; I++) { // 200
         RSCALE[I] = ZERO;
         LSCALE[I] = ZERO;

         WORK[I] = ZERO;
         WORK[I+N] = ZERO;
         WORK[I+2*N] = ZERO;
         WORK[I+3*N] = ZERO;
         WORK[I+4*N] = ZERO;
         WORK[I+5*N] = ZERO;
      } // 200

      // Compute right side vector in resulting linear equations

      BASL = LOG10( SCLFAC );
      for (I = ILO; I <= IHI; I++) { // 240
         for (J = ILO; J <= IHI; J++) { // 230
            TB = B( I, J );
            TA = A( I, J );
            if (TA == ZERO) GO TO 210;
            TA = LOG10( ( TA ).abs() ) / BASL;
            } // 210
            if (TB == ZERO) GO TO 220;
            TB = LOG10( ( TB ).abs() ) / BASL;
            } // 220
            WORK[I+4*N] = WORK( I+4*N ) - TA - TB;
            WORK[J+5*N] = WORK( J+5*N ) - TA - TB;
         } // 230
      } // 240

      COEF = ONE / REAL( 2*NR );
      COEF2 = COEF*COEF;
      COEF5 = HALF*COEF2;
      NRP2 = NR + 2;
      BETA = ZERO;
      IT = 1;

      // Start generalized conjugate gradient iteration

      } // 250

      GAMMA = SDOT( NR, WORK( ILO+4*N ), 1, WORK( ILO+4*N ), 1 ) + SDOT( NR, WORK( ILO+5*N ), 1, WORK( ILO+5*N ), 1 );

      EW = ZERO;
      EWC = ZERO;
      for (I = ILO; I <= IHI; I++) { // 260
         EW = EW + WORK( I+4*N );
         EWC = EWC + WORK( I+5*N );
      } // 260

      GAMMA = COEF*GAMMA - COEF2*( EW**2+EWC**2 ) - COEF5*( EW-EWC )**2;
      if (GAMMA == ZERO) GO TO 350;
      IF( IT != 1 ) BETA = GAMMA / PGAMMA;
      T = COEF5*( EWC-THREE*EW );
      TC = COEF5*( EW-THREE*EWC );

      sscal(NR, BETA, WORK( ILO ), 1 );
      sscal(NR, BETA, WORK( ILO+N ), 1 );

      saxpy(NR, COEF, WORK( ILO+4*N ), 1, WORK( ILO+N ), 1 );
      saxpy(NR, COEF, WORK( ILO+5*N ), 1, WORK( ILO ), 1 );

      for (I = ILO; I <= IHI; I++) { // 270
         WORK[I] = WORK( I ) + TC;
         WORK[I+N] = WORK( I+N ) + T;
      } // 270

      // Apply matrix to vector

      for (I = ILO; I <= IHI; I++) { // 300
         KOUNT = 0;
         SUM = ZERO;
         for (J = ILO; J <= IHI; J++) { // 290
            if( A( I, J ) == ZERO ) GO TO 280;
            KOUNT = KOUNT + 1;
            SUM = SUM + WORK( J );
            } // 280
            if( B( I, J ) == ZERO ) GO TO 290;
            KOUNT = KOUNT + 1;
            SUM = SUM + WORK( J );
         } // 290
         WORK[I+2*N] = double( KOUNT )*WORK( I+N ) + SUM;
      } // 300

      for (J = ILO; J <= IHI; J++) { // 330
         KOUNT = 0;
         SUM = ZERO;
         for (I = ILO; I <= IHI; I++) { // 320
            if( A( I, J ) == ZERO ) GO TO 310;
            KOUNT = KOUNT + 1;
            SUM = SUM + WORK( I+N );
            } // 310
            if( B( I, J ) == ZERO ) GO TO 320;
            KOUNT = KOUNT + 1;
            SUM = SUM + WORK( I+N );
         } // 320
         WORK[J+3*N] = double( KOUNT )*WORK( J ) + SUM;
      } // 330

      SUM = SDOT( NR, WORK( ILO+N ), 1, WORK( ILO+2*N ), 1 ) + SDOT( NR, WORK( ILO ), 1, WORK( ILO+3*N ), 1 );
      ALPHA = GAMMA / SUM;

      // Determine correction to current iteration

      CMAX = ZERO;
      for (I = ILO; I <= IHI; I++) { // 340
         COR = ALPHA*WORK( I+N );
         if( ( COR ).abs() > CMAX ) CMAX = ( COR ).abs();
         LSCALE[I] = LSCALE( I ) + COR;
         COR = ALPHA*WORK( I );
         if( ( COR ).abs() > CMAX ) CMAX = ( COR ).abs();
         RSCALE[I] = RSCALE( I ) + COR;
      } // 340
      if (CMAX < HALF) GO TO 350;

      saxpy(NR, -ALPHA, WORK( ILO+2*N ), 1, WORK( ILO+4*N ), 1 );
      saxpy(NR, -ALPHA, WORK( ILO+3*N ), 1, WORK( ILO+5*N ), 1 );

      PGAMMA = GAMMA;
      IT = IT + 1;
      if (IT <= NRP2) GO TO 250;

      // End generalized conjugate gradient iteration

      } // 350
      SFMIN = SLAMCH( 'S' );
      SFMAX = ONE / SFMIN;
      LSFMIN = INT( LOG10( SFMIN ) / BASL+ONE );
      LSFMAX = INT( LOG10( SFMAX ) / BASL );
      for (I = ILO; I <= IHI; I++) { // 360
         IRAB = ISAMAX( N-ILO+1, A( I, ILO ), LDA );
         RAB = ( A( I, IRAB+ILO-1 ) ).abs();
         IRAB = ISAMAX( N-ILO+1, B( I, ILO ), LDB );
         RAB = max( RAB, ( B( I, IRAB+ILO-1 ) ).abs() );
         LRAB = INT( LOG10( RAB+SFMIN ) / BASL+ONE );
         IR = INT( LSCALE( I ) + sign( HALF, LSCALE( I ) ) );
         IR = min( max( IR, LSFMIN ), LSFMAX, LSFMAX-LRAB );
         LSCALE[I] = SCLFAC**IR;
         ICAB = ISAMAX( IHI, A( 1, I ), 1 );
         CAB = ( A( ICAB, I ) ).abs();
         ICAB = ISAMAX( IHI, B( 1, I ), 1 );
         CAB = max( CAB, ( B( ICAB, I ) ).abs() );
         LCAB = INT( LOG10( CAB+SFMIN ) / BASL+ONE );
         JC = INT( RSCALE( I ) + sign( HALF, RSCALE( I ) ) );
         JC = min( max( JC, LSFMIN ), LSFMAX, LSFMAX-LCAB );
         RSCALE[I] = SCLFAC**JC;
      } // 360

      // Row scaling of matrices A and B

      for (I = ILO; I <= IHI; I++) { // 370
         sscal(N-ILO+1, LSCALE( I ), A( I, ILO ), LDA );
         sscal(N-ILO+1, LSCALE( I ), B( I, ILO ), LDB );
      } // 370

      // Column scaling of matrices A and B

      for (J = ILO; J <= IHI; J++) { // 380
         sscal(IHI, RSCALE( J ), A( 1, J ), 1 );
         sscal(IHI, RSCALE( J ), B( 1, J ), 1 );
      } // 380

      return;
      }
