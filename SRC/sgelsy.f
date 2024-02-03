      void sgelsy(M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS, RANK;
      REAL               RCOND;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      REAL               A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                IMAX, IMIN;
      const              IMAX = 1, IMIN = 2 ;
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IASCL, IBSCL, ISMAX, ISMIN, J, LWKMIN, LWKOPT, MN, NB, NB1, NB2, NB3, NB4;
      REAL               ANRM, BIGNUM, BNRM, C1, C2, S1, S2, SMAX, SMAXPR, SMIN, SMINPR, SMLNUM, WSIZE;
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- REAL               SLAMCH, SLANGE, SROUNDUP_LWORK;
      // EXTERNAL ILAENV, SLAMCH, SLANGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEQP3, SLAIC1, SLASCL, SLASET, SORMQR, SORMRZ, STRSM, STZRZF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      MN = min( M, N );
      ISMIN = MN + 1;
      ISMAX = 2*MN + 1;

      // Test the input arguments.

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, M, N ) ) {
         INFO = -7;
      }

      // Figure out optimal block size

      if ( INFO == 0 ) {
         if ( MN == 0 || NRHS == 0 ) {
            LWKMIN = 1;
            LWKOPT = 1;
         } else {
            NB1 = ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 );
            NB2 = ILAENV( 1, 'SGERQF', ' ', M, N, -1, -1 );
            NB3 = ILAENV( 1, 'SORMQR', ' ', M, N, NRHS, -1 );
            NB4 = ILAENV( 1, 'SORMRQ', ' ', M, N, NRHS, -1 );
            NB = max( NB1, NB2, NB3, NB4 );
            LWKMIN = MN + max( 2*MN, N + 1, MN + NRHS );
            LWKOPT = max( LWKMIN, MN + 2*N + NB*( N + 1 ), 2*MN + NB*NRHS );
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT);

         if ( LWORK < LWKMIN && !LQUERY ) {
            INFO = -12;
         }
      }

      if ( INFO != 0 ) {
         xerbla('SGELSY', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( MN == 0 || NRHS == 0 ) {
         RANK = 0;
         return;
      }

      // Get machine parameters

      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' );
      BIGNUM = ONE / SMLNUM;

      // Scale A, B if max entries outside range [SMLNUM,BIGNUM]

      ANRM = SLANGE( 'M', M, N, A, LDA, WORK );
      IASCL = 0;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         slascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
         IASCL = 1;
      } else if ( ANRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         slascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
         IASCL = 2;
      } else if ( ANRM == ZERO ) {

         // Matrix all zero. Return zero solution.

         slaset('F', max( M, N ), NRHS, ZERO, ZERO, B, LDB );
         RANK = 0;
         GO TO 70;
      }

      BNRM = SLANGE( 'M', M, NRHS, B, LDB, WORK );
      IBSCL = 0;
      if ( BNRM > ZERO && BNRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         slascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 1;
      } else if ( BNRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         slascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 2;
      }

      // Compute QR factorization with column pivoting of A:
         // A * P = Q * R

      sgeqp3(M, N, A, LDA, JPVT, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO );
      WSIZE = MN + WORK( MN+1 );

      // workspace: MN+2*N+NB*(N+1).
      // Details of Householder rotations stored in WORK(1:MN).

      // Determine RANK using incremental condition estimation

      WORK( ISMIN ) = ONE;
      WORK( ISMAX ) = ONE;
      SMAX = ( A( 1, 1 ) ).abs();
      SMIN = SMAX;
      if ( ( A( 1, 1 ) ).abs() == ZERO ) {
         RANK = 0;
         slaset('F', max( M, N ), NRHS, ZERO, ZERO, B, LDB );
         GO TO 70;
      } else {
         RANK = 1;
      }

      } // 10
      if ( RANK < MN ) {
         I = RANK + 1;
         slaic1(IMIN, RANK, WORK( ISMIN ), SMIN, A( 1, I ), A( I, I ), SMINPR, S1, C1 );
         slaic1(IMAX, RANK, WORK( ISMAX ), SMAX, A( 1, I ), A( I, I ), SMAXPR, S2, C2 );

         if ( SMAXPR*RCOND <= SMINPR ) {
            for (I = 1; I <= RANK; I++) { // 20
               WORK( ISMIN+I-1 ) = S1*WORK( ISMIN+I-1 );
               WORK( ISMAX+I-1 ) = S2*WORK( ISMAX+I-1 );
            } // 20
            WORK( ISMIN+RANK ) = C1;
            WORK( ISMAX+RANK ) = C2;
            SMIN = SMINPR;
            SMAX = SMAXPR;
            RANK = RANK + 1;
            GO TO 10;
         }
      }

      // workspace: 3*MN.

      // Logically partition R = [ R11 R12 ]
                              // [  0  R22 ]
      // where R11 = R(1:RANK,1:RANK)

      // [R11,R12] = [ T11, 0 ] * Y

      if (RANK < N) stzrzf( RANK, N, A, LDA, WORK( MN+1 ), WORK( 2*MN+1 ), LWORK-2*MN, INFO );

      // workspace: 2*MN.
      // Details of Householder rotations stored in WORK(MN+1:2*MN)

      // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)

      sormqr('Left', 'Transpose', M, NRHS, MN, A, LDA, WORK( 1 ), B, LDB, WORK( 2*MN+1 ), LWORK-2*MN, INFO );
      WSIZE = max( WSIZE, 2*MN+WORK( 2*MN+1 ) );

      // workspace: 2*MN+NB*NRHS.

      // B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)

      strsm('Left', 'Upper', 'No transpose', 'Non-unit', RANK, NRHS, ONE, A, LDA, B, LDB );

      for (J = 1; J <= NRHS; J++) { // 40
         for (I = RANK + 1; I <= N; I++) { // 30
            B( I, J ) = ZERO;
         } // 30
      } // 40

      // B(1:N,1:NRHS) := Y**T * B(1:N,1:NRHS)

      if ( RANK < N ) {
         sormrz('Left', 'Transpose', N, NRHS, RANK, N-RANK, A, LDA, WORK( MN+1 ), B, LDB, WORK( 2*MN+1 ), LWORK-2*MN, INFO );
      }

      // workspace: 2*MN+NRHS.

      // B(1:N,1:NRHS) := P * B(1:N,1:NRHS)

      for (J = 1; J <= NRHS; J++) { // 60
         for (I = 1; I <= N; I++) { // 50
            WORK( JPVT( I ) ) = B( I, J );
         } // 50
         scopy(N, WORK( 1 ), 1, B( 1, J ), 1 );
      } // 60

      // workspace: N.

      // Undo scaling

      if ( IASCL == 1 ) {
         slascl('G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO );
         slascl('U', 0, 0, SMLNUM, ANRM, RANK, RANK, A, LDA, INFO );
      } else if ( IASCL == 2 ) {
         slascl('G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO );
         slascl('U', 0, 0, BIGNUM, ANRM, RANK, RANK, A, LDA, INFO );
      }
      if ( IBSCL == 1 ) {
         slascl('G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO );
      } else if ( IBSCL == 2 ) {
         slascl('G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO );
      }

      } // 70
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT);

      return;
      }
