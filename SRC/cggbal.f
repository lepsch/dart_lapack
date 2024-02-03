      SUBROUTINE CGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOB;
      int                IHI, ILO, INFO, LDA, LDB, N;
      // ..
      // .. Array Arguments ..
      REAL               LSCALE( * ), RSCALE( * ), WORK( * )
      COMPLEX            A( LDA, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, HALF, ONE
      const              ZERO = 0.0E+0, HALF = 0.5E+0, ONE = 1.0E+0 ;
      REAL               THREE, SCLFAC
      const              THREE = 3.0E+0, SCLFAC = 1.0E+1 ;
      COMPLEX            CZERO
      const              CZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, ICAB, IFLOW, IP1, IR, IRAB, IT, J, JC, JP1, K, KOUNT, L, LCAB, LM1, LRAB, LSFMAX, LSFMIN, M, NR, NRP2;
      REAL               ALPHA, BASL, BETA, CAB, CMAX, COEF, COEF2, COEF5, COR, EW, EWC, GAMMA, PGAMMA, RAB, SFMAX, SFMIN, SUM, T, TA, TB, TC;
      COMPLEX            CDUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ICAMAX;
      REAL               SDOT, SLAMCH
      // EXTERNAL LSAME, ICAMAX, SDOT, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSSCAL, CSWAP, SAXPY, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, INT, LOG10, MAX, MIN, REAL, SIGN
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      if ( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('CGGBAL', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
         ILO = 1
         IHI = N
         RETURN
      }

      if ( N.EQ.1 ) {
         ILO = 1
         IHI = N
         LSCALE( 1 ) = ONE
         RSCALE( 1 ) = ONE
         RETURN
      }

      if ( LSAME( JOB, 'N' ) ) {
         ILO = 1
         IHI = N
         for (I = 1; I <= N; I++) { // 10
            LSCALE( I ) = ONE
            RSCALE( I ) = ONE
   10    CONTINUE
         RETURN
      }

      K = 1
      L = N
      IF( LSAME( JOB, 'S' ) ) GO TO 190

      GO TO 30

      // Permute the matrices A and B to isolate the eigenvalues.

      // Find row with one nonzero in columns 1 through L

   20 CONTINUE
      L = LM1
      IF( L.NE.1 ) GO TO 30

      RSCALE( 1 ) = ONE
      LSCALE( 1 ) = ONE
      GO TO 190

   30 CONTINUE
      LM1 = L - 1
      DO 80 I = L, 1, -1
         for (J = 1; J <= LM1; J++) { // 40
            JP1 = J + 1
            IF( A( I, J ).NE.CZERO .OR. B( I, J ).NE.CZERO ) GO TO 50
   40    CONTINUE
         J = L
         GO TO 70

   50    CONTINUE
         for (J = JP1; J <= L; J++) { // 60
            IF( A( I, J ).NE.CZERO .OR. B( I, J ).NE.CZERO ) GO TO 80
   60    CONTINUE
         J = JP1 - 1

   70    CONTINUE
         M = L
         IFLOW = 1
         GO TO 160
   80 CONTINUE
      GO TO 100

      // Find column with one nonzero in rows K through N

   90 CONTINUE
      K = K + 1

  100 CONTINUE
      for (J = K; J <= L; J++) { // 150
         for (I = K; I <= LM1; I++) { // 110
            IP1 = I + 1
            IF( A( I, J ).NE.CZERO .OR. B( I, J ).NE.CZERO ) GO TO 120
  110    CONTINUE
         I = L
         GO TO 140
  120    CONTINUE
         for (I = IP1; I <= L; I++) { // 130
            IF( A( I, J ).NE.CZERO .OR. B( I, J ).NE.CZERO ) GO TO 150
  130    CONTINUE
         I = IP1 - 1
  140    CONTINUE
         M = K
         IFLOW = 2
         GO TO 160
  150 CONTINUE
      GO TO 190

      // Permute rows M and I

  160 CONTINUE
      LSCALE( M ) = I
      IF( I.EQ.M ) GO TO 170
      cswap(N-K+1, A( I, K ), LDA, A( M, K ), LDA );
      cswap(N-K+1, B( I, K ), LDB, B( M, K ), LDB );

      // Permute columns M and J

  170 CONTINUE
      RSCALE( M ) = J
      IF( J.EQ.M ) GO TO 180
      cswap(L, A( 1, J ), 1, A( 1, M ), 1 );
      cswap(L, B( 1, J ), 1, B( 1, M ), 1 );

  180 CONTINUE
      GO TO ( 20, 90 )IFLOW

  190 CONTINUE
      ILO = K
      IHI = L

      if ( LSAME( JOB, 'P' ) ) {
         for (I = ILO; I <= IHI; I++) { // 195
            LSCALE( I ) = ONE
            RSCALE( I ) = ONE
  195    CONTINUE
         RETURN
      }

      IF( ILO.EQ.IHI ) RETURN

      // Balance the submatrix in rows ILO to IHI.

      NR = IHI - ILO + 1
      for (I = ILO; I <= IHI; I++) { // 200
         RSCALE( I ) = ZERO
         LSCALE( I ) = ZERO

         WORK( I ) = ZERO
         WORK( I+N ) = ZERO
         WORK( I+2*N ) = ZERO
         WORK( I+3*N ) = ZERO
         WORK( I+4*N ) = ZERO
         WORK( I+5*N ) = ZERO
  200 CONTINUE

      // Compute right side vector in resulting linear equations

      BASL = LOG10( SCLFAC )
      for (I = ILO; I <= IHI; I++) { // 240
         for (J = ILO; J <= IHI; J++) { // 230
            if ( A( I, J ).EQ.CZERO ) {
               TA = ZERO
               GO TO 210
            }
            TA = LOG10( CABS1( A( I, J ) ) ) / BASL

  210       CONTINUE
            if ( B( I, J ).EQ.CZERO ) {
               TB = ZERO
               GO TO 220
            }
            TB = LOG10( CABS1( B( I, J ) ) ) / BASL

  220       CONTINUE
            WORK( I+4*N ) = WORK( I+4*N ) - TA - TB
            WORK( J+5*N ) = WORK( J+5*N ) - TA - TB
  230    CONTINUE
  240 CONTINUE

      COEF = ONE / REAL( 2*NR )
      COEF2 = COEF*COEF
      COEF5 = HALF*COEF2
      NRP2 = NR + 2
      BETA = ZERO
      IT = 1

      // Start generalized conjugate gradient iteration

  250 CONTINUE

      GAMMA = SDOT( NR, WORK( ILO+4*N ), 1, WORK( ILO+4*N ), 1 ) + SDOT( NR, WORK( ILO+5*N ), 1, WORK( ILO+5*N ), 1 )

      EW = ZERO
      EWC = ZERO
      for (I = ILO; I <= IHI; I++) { // 260
         EW = EW + WORK( I+4*N )
         EWC = EWC + WORK( I+5*N )
  260 CONTINUE

      GAMMA = COEF*GAMMA - COEF2*( EW**2+EWC**2 ) - COEF5*( EW-EWC )**2
      IF( GAMMA.EQ.ZERO ) GO TO 350       IF( IT.NE.1 ) BETA = GAMMA / PGAMMA
      T = COEF5*( EWC-THREE*EW )
      TC = COEF5*( EW-THREE*EWC )

      sscal(NR, BETA, WORK( ILO ), 1 );
      sscal(NR, BETA, WORK( ILO+N ), 1 );

      saxpy(NR, COEF, WORK( ILO+4*N ), 1, WORK( ILO+N ), 1 );
      saxpy(NR, COEF, WORK( ILO+5*N ), 1, WORK( ILO ), 1 );

      for (I = ILO; I <= IHI; I++) { // 270
         WORK( I ) = WORK( I ) + TC
         WORK( I+N ) = WORK( I+N ) + T
  270 CONTINUE

      // Apply matrix to vector

      for (I = ILO; I <= IHI; I++) { // 300
         KOUNT = 0
         SUM = ZERO
         for (J = ILO; J <= IHI; J++) { // 290
            IF( A( I, J ).EQ.CZERO ) GO TO 280
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( J )
  280       CONTINUE
            IF( B( I, J ).EQ.CZERO ) GO TO 290
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( J )
  290    CONTINUE
         WORK( I+2*N ) = REAL( KOUNT )*WORK( I+N ) + SUM
  300 CONTINUE

      for (J = ILO; J <= IHI; J++) { // 330
         KOUNT = 0
         SUM = ZERO
         for (I = ILO; I <= IHI; I++) { // 320
            IF( A( I, J ).EQ.CZERO ) GO TO 310
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( I+N )
  310       CONTINUE
            IF( B( I, J ).EQ.CZERO ) GO TO 320
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( I+N )
  320    CONTINUE
         WORK( J+3*N ) = REAL( KOUNT )*WORK( J ) + SUM
  330 CONTINUE

      SUM = SDOT( NR, WORK( ILO+N ), 1, WORK( ILO+2*N ), 1 ) + SDOT( NR, WORK( ILO ), 1, WORK( ILO+3*N ), 1 )
      ALPHA = GAMMA / SUM

      // Determine correction to current iteration

      CMAX = ZERO
      for (I = ILO; I <= IHI; I++) { // 340
         COR = ALPHA*WORK( I+N )
         IF( ABS( COR ).GT.CMAX ) CMAX = ABS( COR )
         LSCALE( I ) = LSCALE( I ) + COR
         COR = ALPHA*WORK( I )
         IF( ABS( COR ).GT.CMAX ) CMAX = ABS( COR )
         RSCALE( I ) = RSCALE( I ) + COR
  340 CONTINUE
      IF( CMAX.LT.HALF ) GO TO 350

      saxpy(NR, -ALPHA, WORK( ILO+2*N ), 1, WORK( ILO+4*N ), 1 );
      saxpy(NR, -ALPHA, WORK( ILO+3*N ), 1, WORK( ILO+5*N ), 1 );

      PGAMMA = GAMMA
      IT = IT + 1
      IF( IT.LE.NRP2 ) GO TO 250

      // End generalized conjugate gradient iteration

  350 CONTINUE
      SFMIN = SLAMCH( 'S' )
      SFMAX = ONE / SFMIN
      LSFMIN = INT( LOG10( SFMIN ) / BASL+ONE )
      LSFMAX = INT( LOG10( SFMAX ) / BASL )
      for (I = ILO; I <= IHI; I++) { // 360
         IRAB = ICAMAX( N-ILO+1, A( I, ILO ), LDA )
         RAB = ABS( A( I, IRAB+ILO-1 ) )
         IRAB = ICAMAX( N-ILO+1, B( I, ILO ), LDB )
         RAB = MAX( RAB, ABS( B( I, IRAB+ILO-1 ) ) )
         LRAB = INT( LOG10( RAB+SFMIN ) / BASL+ONE )
         IR = INT( LSCALE( I ) + SIGN( HALF, LSCALE( I ) ) )
         IR = MIN( MAX( IR, LSFMIN ), LSFMAX, LSFMAX-LRAB )
         LSCALE( I ) = SCLFAC**IR
         ICAB = ICAMAX( IHI, A( 1, I ), 1 )
         CAB = ABS( A( ICAB, I ) )
         ICAB = ICAMAX( IHI, B( 1, I ), 1 )
         CAB = MAX( CAB, ABS( B( ICAB, I ) ) )
         LCAB = INT( LOG10( CAB+SFMIN ) / BASL+ONE )
         JC = INT( RSCALE( I ) + SIGN( HALF, RSCALE( I ) ) )
         JC = MIN( MAX( JC, LSFMIN ), LSFMAX, LSFMAX-LCAB )
         RSCALE( I ) = SCLFAC**JC
  360 CONTINUE

      // Row scaling of matrices A and B

      for (I = ILO; I <= IHI; I++) { // 370
         csscal(N-ILO+1, LSCALE( I ), A( I, ILO ), LDA );
         csscal(N-ILO+1, LSCALE( I ), B( I, ILO ), LDB );
  370 CONTINUE

      // Column scaling of matrices A and B

      for (J = ILO; J <= IHI; J++) { // 380
         csscal(IHI, RSCALE( J ), A( 1, J ), 1 );
         csscal(IHI, RSCALE( J ), B( 1, J ), 1 );
  380 CONTINUE

      RETURN

      // End of CGGBAL

      }
