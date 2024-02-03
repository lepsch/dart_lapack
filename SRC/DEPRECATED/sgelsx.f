      SUBROUTINE SGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, M, N, NRHS, RANK;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      REAL               A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                IMAX, IMIN;
      const              IMAX = 1, IMIN = 2 ;
      REAL               ZERO, ONE, DONE, NTDONE
      const              ZERO = 0.0E0, ONE = 1.0E0, DONE = ZERO, NTDONE = ONE ;
      // ..
      // .. Local Scalars ..
      int                I, IASCL, IBSCL, ISMAX, ISMIN, J, K, MN;
      REAL               ANRM, BIGNUM, BNRM, C1, C2, S1, S2, SMAX, SMAXPR, SMIN, SMINPR, SMLNUM, T1, T2
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEQPF, SLAIC1, SLASCL, SLASET, SLATZM, SORM2R, STRSM, STZRQF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      MN = MIN( M, N )
      ISMIN = MN + 1
      ISMAX = 2*MN + 1

      // Test the input arguments.

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, M, N ) ) {
         INFO = -7
      }

      if ( INFO.NE.0 ) {
         xerbla('SGELSX', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( MIN( M, N, NRHS ).EQ.0 ) {
         RANK = 0
         RETURN
      }

      // Get machine parameters

      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale A, B if max elements outside range [SMLNUM,BIGNUM]

      ANRM = SLANGE( 'M', M, N, A, LDA, WORK )
      IASCL = 0
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         slascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
         IASCL = 1
      } else if ( ANRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         slascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
         IASCL = 2
      } else if ( ANRM.EQ.ZERO ) {

         // Matrix all zero. Return zero solution.

         slaset('F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB );
         RANK = 0
         GO TO 100
      }

      BNRM = SLANGE( 'M', M, NRHS, B, LDB, WORK )
      IBSCL = 0
      if ( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         slascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 1
      } else if ( BNRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         slascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 2
      }

      // Compute QR factorization with column pivoting of A:
         // A * P = Q * R

      sgeqpf(M, N, A, LDA, JPVT, WORK( 1 ), WORK( MN+1 ), INFO );

      // workspace 3*N. Details of Householder rotations stored
      // in WORK(1:MN).

      // Determine RANK using incremental condition estimation

      WORK( ISMIN ) = ONE
      WORK( ISMAX ) = ONE
      SMAX = ABS( A( 1, 1 ) )
      SMIN = SMAX
      if ( ABS( A( 1, 1 ) ).EQ.ZERO ) {
         RANK = 0
         slaset('F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB );
         GO TO 100
      } else {
         RANK = 1
      }

      } // 10
      if ( RANK.LT.MN ) {
         I = RANK + 1
         slaic1(IMIN, RANK, WORK( ISMIN ), SMIN, A( 1, I ), A( I, I ), SMINPR, S1, C1 )          CALL SLAIC1( IMAX, RANK, WORK( ISMAX ), SMAX, A( 1, I ), A( I, I ), SMAXPR, S2, C2 );

         if ( SMAXPR*RCOND.LE.SMINPR ) {
            for (I = 1; I <= RANK; I++) { // 20
               WORK( ISMIN+I-1 ) = S1*WORK( ISMIN+I-1 )
               WORK( ISMAX+I-1 ) = S2*WORK( ISMAX+I-1 )
            } // 20
            WORK( ISMIN+RANK ) = C1
            WORK( ISMAX+RANK ) = C2
            SMIN = SMINPR
            SMAX = SMAXPR
            RANK = RANK + 1
            GO TO 10
         }
      }

      // Logically partition R = [ R11 R12 ]
                              // [  0  R22 ]
      // where R11 = R(1:RANK,1:RANK)

      // [R11,R12] = [ T11, 0 ] * Y

      IF( RANK.LT.N ) CALL STZRQF( RANK, N, A, LDA, WORK( MN+1 ), INFO )

      // Details of Householder rotations stored in WORK(MN+1:2*MN)

      // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)

      sorm2r('Left', 'Transpose', M, NRHS, MN, A, LDA, WORK( 1 ), B, LDB, WORK( 2*MN+1 ), INFO );

      // workspace NRHS

      // B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)

      strsm('Left', 'Upper', 'No transpose', 'Non-unit', RANK, NRHS, ONE, A, LDA, B, LDB );

      DO 40 I = RANK + 1, N
         for (J = 1; J <= NRHS; J++) { // 30
            B( I, J ) = ZERO
         } // 30
      } // 40

      // B(1:N,1:NRHS) := Y**T * B(1:N,1:NRHS)

      if ( RANK.LT.N ) {
         for (I = 1; I <= RANK; I++) { // 50
            slatzm('Left', N-RANK+1, NRHS, A( I, RANK+1 ), LDA, WORK( MN+I ), B( I, 1 ), B( RANK+1, 1 ), LDB, WORK( 2*MN+1 ) );
         } // 50
      }

      // workspace NRHS

      // B(1:N,1:NRHS) := P * B(1:N,1:NRHS)

      for (J = 1; J <= NRHS; J++) { // 90
         for (I = 1; I <= N; I++) { // 60
            WORK( 2*MN+I ) = NTDONE
         } // 60
         for (I = 1; I <= N; I++) { // 80
            if ( WORK( 2*MN+I ).EQ.NTDONE ) {
               if ( JPVT( I ).NE.I ) {
                  K = I
                  T1 = B( K, J )
                  T2 = B( JPVT( K ), J )
                  } // 70
                  B( JPVT( K ), J ) = T1
                  WORK( 2*MN+K ) = DONE
                  T1 = T2
                  K = JPVT( K )
                  T2 = B( JPVT( K ), J )
                  IF( JPVT( K ).NE.I ) GO TO 70
                  B( I, J ) = T1
                  WORK( 2*MN+K ) = DONE
               }
            }
         } // 80
      } // 90

      // Undo scaling

      if ( IASCL.EQ.1 ) {
         slascl('G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO );
         slascl('U', 0, 0, SMLNUM, ANRM, RANK, RANK, A, LDA, INFO );
      } else if ( IASCL.EQ.2 ) {
         slascl('G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO );
         slascl('U', 0, 0, BIGNUM, ANRM, RANK, RANK, A, LDA, INFO );
      }
      if ( IBSCL.EQ.1 ) {
         slascl('G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO );
      } else if ( IBSCL.EQ.2 ) {
         slascl('G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO );
      }

      } // 100

      RETURN

      // End of SGELSX

      }
