      SUBROUTINE DTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, WR, WI, M, S, SEP, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPQ, JOB;
      int                INFO, LDQ, LDT, LIWORK, LWORK, M, N;
      double             S, SEP;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      int                IWORK( * );
      double             Q( LDQ, * ), T( LDT, * ), WI( * ), WORK( * ), WR( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, PAIR, SWAP, WANTBH, WANTQ, WANTS, WANTSP;
      int                IERR, K, KASE, KK, KS, LIWMIN, LWMIN, N1, N2, NN;
      double             EST, RNORM, SCALE;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLANGE;
      // EXTERNAL LSAME, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, DLACPY, DTREXC, DTRSYL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
      WANTQ = LSAME( COMPQ, 'V' )

      INFO = 0
      LQUERY = ( LWORK == -1 )
      if ( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.WANTS .AND. .NOT.WANTSP ) {
         INFO = -1
      } else if ( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.WANTQ ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDT.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) {
         INFO = -8
      } else {

         // Set M to the dimension of the specified invariant subspace,
         // and test LWORK and LIWORK.

         M = 0
         PAIR = false;
         for (K = 1; K <= N; K++) { // 10
            if ( PAIR ) {
               PAIR = false;
            } else {
               if ( K.LT.N ) {
                  if ( T( K+1, K ) == ZERO ) {
                     IF( SELECT( K ) ) M = M + 1
                  } else {
                     PAIR = true;
                     IF( SELECT( K ) .OR. SELECT( K+1 ) ) M = M + 2
                  }
               } else {
                  IF( SELECT( N ) ) M = M + 1
               }
            }
         } // 10

         N1 = M
         N2 = N - M
         NN = N1*N2

         if ( WANTSP ) {
            LWMIN = MAX( 1, 2*NN )
            LIWMIN = MAX( 1, NN )
         } else if ( LSAME( JOB, 'N' ) ) {
            LWMIN = MAX( 1, N )
            LIWMIN = 1
         } else if ( LSAME( JOB, 'E' ) ) {
            LWMIN = MAX( 1, NN )
            LIWMIN = 1
         }

         if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -15
         } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
            INFO = -17
         }
      }

      if ( INFO == 0 ) {
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      }

      if ( INFO != 0 ) {
         xerbla('DTRSEN', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible.

      if ( M == N .OR. M == 0 ) {
         if (WANTS) S = ONE          IF( WANTSP ) SEP = DLANGE( '1', N, N, T, LDT, WORK );
         GO TO 40
      }

      // Collect the selected blocks at the top-left corner of T.

      KS = 0
      PAIR = false;
      for (K = 1; K <= N; K++) { // 20
         if ( PAIR ) {
            PAIR = false;
         } else {
            SWAP = SELECT( K )
            if ( K.LT.N ) {
               if ( T( K+1, K ) != ZERO ) {
                  PAIR = true;
                  SWAP = SWAP .OR. SELECT( K+1 )
               }
            }
            if ( SWAP ) {
               KS = KS + 1

               // Swap the K-th block to position KS.

               IERR = 0
               KK = K
               if (K != KS) CALL DTREXC( COMPQ, N, T, LDT, Q, LDQ, KK, KS, WORK, IERR );
               if ( IERR == 1 .OR. IERR == 2 ) {

                  // Blocks too close to swap: exit.

                  INFO = 1
                  if (WANTS) S = ZERO                   IF( WANTSP ) SEP = ZERO;
                  GO TO 40
               }
               if (PAIR) KS = KS + 1;
            }
         }
      } // 20

      if ( WANTS ) {

         // Solve Sylvester equation for R:

            // T11*R - R*T22 = scale*T12

         dlacpy('F', N1, N2, T( 1, N1+1 ), LDT, WORK, N1 );
         dtrsyl('N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR );

         // Estimate the reciprocal of the condition number of the cluster
         // of eigenvalues.

         RNORM = DLANGE( 'F', N1, N2, WORK, N1, WORK )
         if ( RNORM == ZERO ) {
            S = ONE
         } else {
            S = SCALE / ( SQRT( SCALE*SCALE / RNORM+RNORM )* SQRT( RNORM ) )
         }
      }

      if ( WANTSP ) {

         // Estimate sep(T11,T22).

         EST = ZERO
         KASE = 0
         } // 30
         dlacn2(NN, WORK( NN+1 ), WORK, IWORK, EST, KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == 1 ) {

               // Solve  T11*R - R*T22 = scale*X.

               dtrsyl('N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR );
            } else {

               // Solve T11**T*R - R*T22**T = scale*X.

               dtrsyl('T', 'T', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR );
            }
            GO TO 30
         }

         SEP = SCALE / EST
      }

      } // 40

      // Store the output eigenvalues in WR and WI.

      for (K = 1; K <= N; K++) { // 50
         WR( K ) = T( K, K )
         WI( K ) = ZERO
      } // 50
      for (K = 1; K <= N - 1; K++) { // 60
         if ( T( K+1, K ) != ZERO ) {
            WI( K ) = SQRT( ABS( T( K, K+1 ) ) )* SQRT( ABS( T( K+1, K ) ) )
            WI( K+1 ) = -WI( K )
         }
      } // 60

      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of DTRSEN

      }
