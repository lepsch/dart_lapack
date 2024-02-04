      void dgelsx(M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, M, N, NRHS, RANK;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      double             A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                IMAX, IMIN;
      const              IMAX = 1, IMIN = 2 ;
      double             ZERO, ONE, DONE, NTDONE;
      const              ZERO = 0.0, ONE = 1.0, DONE = ZERO, NTDONE = ONE ;
      // ..
      // .. Local Scalars ..
      int                I, IASCL, IBSCL, ISMAX, ISMIN, J, K, MN;
      double             ANRM, BIGNUM, BNRM, C1, C2, S1, S2, SMAX, SMAXPR, SMIN, SMINPR, SMLNUM, T1, T2;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQPF, DLAIC1, DLASCL, DLASET, DLATZM, DORM2R, DTRSM, DTZRQF, XERBLA
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

      if ( INFO != 0 ) {
         xerbla('DGELSX', -INFO );
         return;
      }

      // Quick return if possible

      if ( min( M, N, NRHS ) == 0 ) {
         RANK = 0;
         return;
      }

      // Get machine parameters

      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' );
      BIGNUM = ONE / SMLNUM;

      // Scale A, B if max elements outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', M, N, A, LDA, WORK );
      IASCL = 0;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         dlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
         IASCL = 1;
      } else if ( ANRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         dlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
         IASCL = 2;
      } else if ( ANRM == ZERO ) {

         // Matrix all zero. Return zero solution.

         dlaset('F', max( M, N ), NRHS, ZERO, ZERO, B, LDB );
         RANK = 0;
         GO TO 100;
      }

      BNRM = DLANGE( 'M', M, NRHS, B, LDB, WORK );
      IBSCL = 0;
      if ( BNRM > ZERO && BNRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         dlascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 1;
      } else if ( BNRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         dlascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 2;
      }

      // Compute QR factorization with column pivoting of A:
         // A * P = Q * R

      dgeqpf(M, N, A, LDA, JPVT, WORK( 1 ), WORK( MN+1 ), INFO );

      // workspace 3*N. Details of Householder rotations stored
      // in WORK(1:MN).

      // Determine RANK using incremental condition estimation

      WORK[ISMIN] = ONE;
      WORK[ISMAX] = ONE;
      SMAX = ( A( 1, 1 ) ).abs();
      SMIN = SMAX;
      if ( ( A( 1, 1 ) ).abs() == ZERO ) {
         RANK = 0;
         dlaset('F', max( M, N ), NRHS, ZERO, ZERO, B, LDB );
         GO TO 100;
      } else {
         RANK = 1;
      }

      } // 10
      if ( RANK < MN ) {
         I = RANK + 1;
         dlaic1(IMIN, RANK, WORK( ISMIN ), SMIN, A( 1, I ), A( I, I ), SMINPR, S1, C1 );
         dlaic1(IMAX, RANK, WORK( ISMAX ), SMAX, A( 1, I ), A( I, I ), SMAXPR, S2, C2 );

         if ( SMAXPR*RCOND <= SMINPR ) {
            for (I = 1; I <= RANK; I++) { // 20
               WORK[ISMIN+I-1] = S1*WORK( ISMIN+I-1 );
               WORK[ISMAX+I-1] = S2*WORK( ISMAX+I-1 );
            } // 20
            WORK[ISMIN+RANK] = C1;
            WORK[ISMAX+RANK] = C2;
            SMIN = SMINPR;
            SMAX = SMAXPR;
            RANK = RANK + 1;
            GO TO 10;
         }
      }

      // Logically partition R = [ R11 R12 ]
                              // [  0  R22 ]
      // where R11 = R(1:RANK,1:RANK)

      // [R11,R12] = [ T11, 0 ] * Y

      if (RANK < N) dtzrqf( RANK, N, A, LDA, WORK( MN+1 ), INFO );

      // Details of Householder rotations stored in WORK(MN+1:2*MN)

      // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)

      dorm2r('Left', 'Transpose', M, NRHS, MN, A, LDA, WORK( 1 ), B, LDB, WORK( 2*MN+1 ), INFO );

      // workspace NRHS

      // B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)

      dtrsm('Left', 'Upper', 'No transpose', 'Non-unit', RANK, NRHS, ONE, A, LDA, B, LDB );

      for (I = RANK + 1; I <= N; I++) { // 40
         for (J = 1; J <= NRHS; J++) { // 30
            B[I, J] = ZERO;
         } // 30
      } // 40

      // B(1:N,1:NRHS) := Y**T * B(1:N,1:NRHS)

      if ( RANK < N ) {
         for (I = 1; I <= RANK; I++) { // 50
            dlatzm('Left', N-RANK+1, NRHS, A( I, RANK+1 ), LDA, WORK( MN+I ), B( I, 1 ), B( RANK+1, 1 ), LDB, WORK( 2*MN+1 ) );
         } // 50
      }

      // workspace NRHS

      // B(1:N,1:NRHS) := P * B(1:N,1:NRHS)

      for (J = 1; J <= NRHS; J++) { // 90
         for (I = 1; I <= N; I++) { // 60
            WORK[2*MN+I] = NTDONE;
         } // 60
         for (I = 1; I <= N; I++) { // 80
            if ( WORK( 2*MN+I ) == NTDONE ) {
               if ( JPVT( I ) != I ) {
                  K = I;
                  T1 = B( K, J );
                  T2 = B( JPVT( K ), J );
                  } // 70
                  B[JPVT( K ), J] = T1;
                  WORK[2*MN+K] = DONE;
                  T1 = T2;
                  K = JPVT( K );
                  T2 = B( JPVT( K ), J );
                  if( JPVT( K ) != I ) GO TO 70;
                  B[I, J] = T1;
                  WORK[2*MN+K] = DONE;
               }
            }
         } // 80
      } // 90

      // Undo scaling

      if ( IASCL == 1 ) {
         dlascl('G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO );
         dlascl('U', 0, 0, SMLNUM, ANRM, RANK, RANK, A, LDA, INFO );
      } else if ( IASCL == 2 ) {
         dlascl('G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO );
         dlascl('U', 0, 0, BIGNUM, ANRM, RANK, RANK, A, LDA, INFO );
      }
      if ( IBSCL == 1 ) {
         dlascl('G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO );
      } else if ( IBSCL == 2 ) {
         dlascl('G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO );
      }

      } // 100

      return;
      }
