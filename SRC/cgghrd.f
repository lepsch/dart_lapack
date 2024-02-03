      SUBROUTINE CGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q, LDQ, Z, LDZ, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPQ, COMPZ;
      int                IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            CONE, CZERO;
      const              CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ILQ, ILZ;
      int                ICOMPQ, ICOMPZ, JCOL, JROW;
      REAL               C;
      COMPLEX            CTEMP, S;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARTG, CLASET, CROT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX
      // ..
      // .. Executable Statements ..

      // Decode COMPQ

      if ( LSAME( COMPQ, 'N' ) ) {
         ILQ = false;
         ICOMPQ = 1;
      } else if ( LSAME( COMPQ, 'V' ) ) {
         ILQ = true;
         ICOMPQ = 2;
      } else if ( LSAME( COMPQ, 'I' ) ) {
         ILQ = true;
         ICOMPQ = 3;
      } else {
         ICOMPQ = 0;
      }

      // Decode COMPZ

      if ( LSAME( COMPZ, 'N' ) ) {
         ILZ = false;
         ICOMPZ = 1;
      } else if ( LSAME( COMPZ, 'V' ) ) {
         ILZ = true;
         ICOMPZ = 2;
      } else if ( LSAME( COMPZ, 'I' ) ) {
         ILZ = true;
         ICOMPZ = 3;
      } else {
         ICOMPZ = 0;
      }

      // Test the input parameters.

      INFO = 0;
      if ( ICOMPQ <= 0 ) {
         INFO = -1;
      } else if ( ICOMPZ <= 0 ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( ILO < 1 ) {
         INFO = -4;
      } else if ( IHI > N || IHI < ILO-1 ) {
         INFO = -5;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -7;
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -9;
      } else if ( ( ILQ && LDQ < N ) || LDQ < 1 ) {
         INFO = -11;
      } else if ( ( ILZ && LDZ < N ) || LDZ < 1 ) {
         INFO = -13;
      }
      if ( INFO != 0 ) {
         xerbla('CGGHRD', -INFO );
         RETURN;
      }

      // Initialize Q and Z if desired.

      if (ICOMPQ == 3) CALL CLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )       IF( ICOMPZ == 3 ) CALL CLASET( 'Full', N, N, CZERO, CONE, Z, LDZ );

      // Quick return if possible

      if (N <= 1) RETURN;

      // Zero out lower triangle of B

      for (JCOL = 1; JCOL <= N - 1; JCOL++) { // 20
         for (JROW = JCOL + 1; JROW <= N; JROW++) { // 10
            B( JROW, JCOL ) = CZERO;
         } // 10
      } // 20

      // Reduce A and B

      for (JCOL = ILO; JCOL <= IHI - 2; JCOL++) { // 40

         DO 30 JROW = IHI, JCOL + 2, -1;

            // Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)

            CTEMP = A( JROW-1, JCOL );
            clartg(CTEMP, A( JROW, JCOL ), C, S, A( JROW-1, JCOL ) );
            A( JROW, JCOL ) = CZERO;
            crot(N-JCOL, A( JROW-1, JCOL+1 ), LDA, A( JROW, JCOL+1 ), LDA, C, S );
            crot(N+2-JROW, B( JROW-1, JROW-1 ), LDB, B( JROW, JROW-1 ), LDB, C, S )             IF( ILQ ) CALL CROT( N, Q( 1, JROW-1 ), 1, Q( 1, JROW ), 1, C, CONJG( S ) );

            // Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)

            CTEMP = B( JROW, JROW );
            clartg(CTEMP, B( JROW, JROW-1 ), C, S, B( JROW, JROW ) );
            B( JROW, JROW-1 ) = CZERO;
            crot(IHI, A( 1, JROW ), 1, A( 1, JROW-1 ), 1, C, S );
            crot(JROW-1, B( 1, JROW ), 1, B( 1, JROW-1 ), 1, C, S )             IF( ILZ ) CALL CROT( N, Z( 1, JROW ), 1, Z( 1, JROW-1 ), 1, C, S );
         } // 30
      } // 40

      RETURN;

      // End of CGGHRD

      }
