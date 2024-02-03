      SUBROUTINE SGEQP3RK( M, N, NRHS, KMAX, ABSTOL, RELTOL, A, LDA, K, MAXC2NRMK, RELMAXC2NRMK, JPIV, TAU, WORK, LWORK, IWORK, INFO )
      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, KF, KMAX, LDA, LWORK, M, N, NRHS;
      REAL               ABSTOL,  MAXC2NRMK, RELMAXC2NRMK, RELTOL
      // ..
      // .. Array Arguments ..
      int                IWORK( * ), JPIV( * );
      REAL               A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                INB, INBMIN, IXOVER;
      PARAMETER          ( INB = 1, INBMIN = 2, IXOVER = 3 )
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               LQUERY, DONE;
      int                IINFO, IOFFSET, IWS, J, JB, JBF, JMAXB, JMAX, JMAXC2NRM, KP1, LWKOPT, MINMN, N_SUB, NB, NBMIN, NX;
      REAL               EPS, HUGEVAL, MAXC2NRM, SAFMIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAQP2RK, SLAQP3RK, XERBLA
      // ..
      // .. External Functions ..
      bool               SISNAN;
      int                ISAMAX, ILAENV;
      REAL               SLAMCH, SNRM2, SROUNDUP_LWORK
      // EXTERNAL SISNAN, SLAMCH, SNRM2, ISAMAX, ILAENV, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test input arguments
      // ====================

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( KMAX.LT.0 ) THEN
         INFO = -4
      ELSE IF( SISNAN( ABSTOL ) ) THEN
         INFO = -5
      ELSE IF( SISNAN( RELTOL ) ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -8
      END IF

      // If the input parameters M, N, NRHS, KMAX, LDA are valid:
        // a) Test the input workspace size LWORK for the minimum
           // size requirement IWS.
        // b) Determine the optimal block size NB and optimal
           // workspace size LWKOPT to be returned in WORK(1)
           // in case of (1) LWORK < IWS, (2) LQUERY = .TRUE.,
           // (3) when routine exits.
      // Here, IWS is the miminum workspace required for unblocked
      // code.

      IF( INFO.EQ.0 ) THEN
         MINMN = MIN( M, N )
         IF( MINMN.EQ.0 ) THEN
            IWS = 1
            LWKOPT = 1
         ELSE

            // Minimal workspace size in case of using only unblocked
            // BLAS 2 code in SLAQP2RK.
            // 1) SGEQP3RK and SLAQP2RK: 2*N to store full and partial
               // column 2-norms.
            // 2) SLAQP2RK: N+NRHS-1 to use in WORK array that is used
               // in SLARF subroutine inside SLAQP2RK to apply an
               // elementary reflector from the left.
            // TOTAL_WORK_SIZE = 3*N + NRHS - 1

            IWS = 3*N + NRHS - 1

            // Assign to NB optimal block size.

            NB = ILAENV( INB, 'SGEQP3RK', ' ', M, N, -1, -1 )

            // A formula for the optimal workspace size in case of using
            // both unblocked BLAS 2 in SLAQP2RK and blocked BLAS 3 code
            // in SLAQP3RK.
            // 1) SGEQP3RK, SLAQP2RK, SLAQP3RK: 2*N to store full and
               // partial column 2-norms.
            // 2) SLAQP2RK: N+NRHS-1 to use in WORK array that is used
               // in SLARF subroutine to apply an elementary reflector
               // from the left.
            // 3) SLAQP3RK: NB*(N+NRHS) to use in the work array F that
               // is used to apply a block reflector from
              t // he left.
            // 4) SLAQP3RK: NB to use in the auxilixary array AUX.
            // Sizes (2) and ((3) + (4)) should intersect, therefore
            // TOTAL_WORK_SIZE = 2*N + NB*( N+NRHS+1 ), given NBMIN=2.

            LWKOPT = 2*N + NB*( N+NRHS+1 )
         END IF
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

         IF( ( LWORK.LT.IWS ) .AND. .NOT.LQUERY ) THEN
            INFO = -15
         END IF
      END IF

       // NOTE: The optimal workspace size is returned in WORK(1), if
            t // he input parameters M, N, NRHS, KMAX, LDA are valid.

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEQP3RK', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible for M=0 or N=0.

      IF( MINMN.EQ.0 ) THEN
         K = 0
         MAXC2NRMK = ZERO
         RELMAXC2NRMK = ZERO
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
         RETURN
      END IF

      // ==================================================================

      // Initialize column pivot array JPIV.

      DO J = 1, N
         JPIV( J ) = J
      END DO

      // ==================================================================

      // Initialize storage for partial and exact column 2-norms.
      // a) The elements WORK(1:N) are used to store partial column
         // 2-norms of the matrix A, and may decrease in each computation
         // step; initialize to the values of complete columns 2-norms.
      // b) The elements WORK(N+1:2*N) are used to store complete column
         // 2-norms of the matrix A, they are not changed during the
         // computation; initialize the values of complete columns 2-norms.

      DO J = 1, N
         WORK( J ) = SNRM2( M, A( 1, J ), 1 )
         WORK( N+J ) = WORK( J )
      END DO

      // ==================================================================

      // Compute the pivot column index and the maximum column 2-norm
      // for the whole original matrix stored in A(1:M,1:N).

      KP1 = ISAMAX( N, WORK( 1 ), 1 )
      MAXC2NRM = WORK( KP1 )

      // ==================================================================.

      IF( SISNAN( MAXC2NRM ) ) THEN

         // Check if the matrix A contains NaN, set INFO parameter
        t // o the column number where the first NaN is found and return
         // from the routine.

         K = 0
         INFO = KP1

         // Set MAXC2NRMK and  RELMAXC2NRMK to NaN.

         MAXC2NRMK = MAXC2NRM
         RELMAXC2NRMK = MAXC2NRM

         // Array TAU is not set and contains undefined elements.

         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
         RETURN
      END IF

      // ===================================================================

      IF( MAXC2NRM.EQ.ZERO ) THEN

         // Check is the matrix A is a zero matrix, set array TAU and
         // return from the routine.

         K = 0
         MAXC2NRMK = ZERO
         RELMAXC2NRMK = ZERO

         DO J = 1, MINMN
            TAU( J ) = ZERO
         END DO

         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
         RETURN

      END IF

      // ===================================================================

      HUGEVAL = SLAMCH( 'Overflow' )

      IF( MAXC2NRM.GT.HUGEVAL ) THEN

         // Check if the matrix A contains +Inf or -Inf, set INFO parameter
        t // o the column number, where the first +/-Inf  is found plus N,
         // and continue the computation.

         INFO = N + KP1

      END IF

      // ==================================================================

      // Quick return if possible for the case when the first
      // stopping criterion is satisfied, i.e. KMAX = 0.

      IF( KMAX.EQ.0 ) THEN
         K = 0
         MAXC2NRMK = MAXC2NRM
         RELMAXC2NRMK = ONE
         DO J = 1, MINMN
            TAU( J ) = ZERO
         END DO
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
         RETURN
      END IF

      // ==================================================================

      EPS = SLAMCH('Epsilon')

      // Adjust ABSTOL

      IF( ABSTOL.GE.ZERO ) THEN
         SAFMIN = SLAMCH('Safe minimum')
         ABSTOL = MAX( ABSTOL, TWO*SAFMIN )
      END IF

      // Adjust RELTOL

      IF( RELTOL.GE.ZERO ) THEN
         RELTOL = MAX( RELTOL, EPS )
      END IF

      // ===================================================================

      // JMAX is the maximum index of the column to be factorized,
      // which is also limited by the first stopping criterion KMAX.

      JMAX = MIN( KMAX, MINMN )

      // ===================================================================

      // Quick return if possible for the case when the second or third
      // stopping criterion for the whole original matrix is satified,
      // i.e. MAXC2NRM <= ABSTOL or RELMAXC2NRM <= RELTOL
      // (which is ONE <= RELTOL).

      IF( MAXC2NRM.LE.ABSTOL .OR. ONE.LE.RELTOL ) THEN

         K = 0
         MAXC2NRMK = MAXC2NRM
         RELMAXC2NRMK = ONE

         DO J = 1, MINMN
            TAU( J ) = ZERO
         END DO

         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
         RETURN
      END IF

      // ==================================================================
      // Factorize columns
      // ==================================================================

      // Determine the block size.

      NBMIN = 2
      NX = 0

      IF( ( NB.GT.1 ) .AND. ( NB.LT.MINMN ) ) THEN

         // Determine when to cross over from blocked to unblocked code.
         // (for N less than NX, unblocked code should be used).

         NX = MAX( 0, ILAENV( IXOVER, 'SGEQP3RK', ' ', M, N, -1, -1 ))

         IF( NX.LT.MINMN ) THEN

            // Determine if workspace is large enough for blocked code.

            IF( LWORK.LT.LWKOPT ) THEN

               // Not enough workspace to use optimal block size that
               // is currently stored in NB.
               // Reduce NB and determine the minimum value of NB.

               NB = ( LWORK-2*N ) / ( N+1 )
               NBMIN = MAX( 2, ILAENV( INBMIN, 'SGEQP3RK', ' ', M, N, -1, -1 ) )

            END IF
         END IF
      END IF

      // ==================================================================

      // DONE is the boolean flag to rerpresent the case when the
      // factorization completed in the block factorization routine,
      // before the end of the block.

      DONE = .FALSE.

      // J is the column index.

      J = 1

      // (1) Use blocked code initially.

      // JMAXB is the maximum column index of the block, when the
      // blocked code is used, is also limited by the first stopping
      // criterion KMAX.

      JMAXB = MIN( KMAX, MINMN - NX )

      IF( NB.GE.NBMIN .AND. NB.LT.JMAX .AND. JMAXB.GT.0 ) THEN

         // Loop over the column blocks of the matrix A(1:M,1:JMAXB). Here:
         // J   is the column index of a column block;
         // JB  is the column block size to pass to block factorization
             // routine in a loop step;
         // JBF is the number of columns that were actually factorized
            t // hat was returned by the block factorization routine
             // in a loop step, JBF <= JB;
         // N_SUB is the number of columns in the submatrix;
         // IOFFSET is the number of rows that should not be factorized.

         DO WHILE( J.LE.JMAXB )

            JB = MIN( NB, JMAXB-J+1 )
            N_SUB = N-J+1
            IOFFSET = J-1

            // Factorize JB columns among the columns A(J:N).

            CALL SLAQP3RK( M, N_SUB, NRHS, IOFFSET, JB, ABSTOL, RELTOL, KP1, MAXC2NRM, A( 1, J ), LDA, DONE, JBF, MAXC2NRMK, RELMAXC2NRMK, JPIV( J ), TAU( J ), WORK( J ), WORK( N+J ), WORK( 2*N+1 ), WORK( 2*N+JB+1 ), N+NRHS-J+1, IWORK, IINFO )

            // Set INFO on the first occurence of Inf.

            IF( IINFO.GT.N_SUB .AND. INFO.EQ.0 ) THEN
               INFO = 2*IOFFSET + IINFO
            END IF

            IF( DONE ) THEN

               // Either the submatrix is zero before the end of the
               // column block, or ABSTOL or RELTOL criterion is
               // satisfied before the end of the column block, we can
               // return from the routine. Perform the following before
               // returning:
                 // a) Set the number of factorized columns K,
                    // K = IOFFSET + JBF from the last call of blocked
                    // routine.
                 // NOTE: 1) MAXC2NRMK and RELMAXC2NRMK are returned
                          // by the block factorization routine;
                       // 2) The remaining TAUs are set to ZERO by the
                          // block factorization routine.

               K = IOFFSET + JBF

               // Set INFO on the first occurrence of NaN, NaN takes
               // prcedence over Inf.

               IF( IINFO.LE.N_SUB .AND. IINFO.GT.0 ) THEN
                  INFO = IOFFSET + IINFO
               END IF

               // Return from the routine.

               WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

               RETURN

            END IF

            J = J + JBF

         END DO

      END IF

      // Use unblocked code to factor the last or only block.
      // J = JMAX+1 means we factorized the maximum possible number of
      // columns, that is in ELSE clause we need to compute
     t // he MAXC2NORM and RELMAXC2NORM to return after we processed
     t // he blocks.

      IF( J.LE.JMAX ) THEN

         // N_SUB is the number of columns in the submatrix;
         // IOFFSET is the number of rows that should not be factorized.

         N_SUB = N-J+1
         IOFFSET = J-1

         CALL SLAQP2RK( M, N_SUB, NRHS, IOFFSET, JMAX-J+1, ABSTOL, RELTOL, KP1, MAXC2NRM, A( 1, J ), LDA, KF, MAXC2NRMK, RELMAXC2NRMK, JPIV( J ), TAU( J ), WORK( J ), WORK( N+J ), WORK( 2*N+1 ), IINFO )

         // ABSTOL or RELTOL criterion is satisfied when the number of
        t // he factorized columns KF is smaller then the  number
         // of columns JMAX-J+1 supplied to be factorized by the
         // unblocked routine, we can return from
        t // he routine. Perform the following before returning:
            // a) Set the number of factorized columns K,
            // b) MAXC2NRMK and RELMAXC2NRMK are returned by the
               // unblocked factorization routine above.

         K = J - 1 + KF

         // Set INFO on the first exception occurence.

         // Set INFO on the first exception occurence of Inf or NaN,
         // (NaN takes precedence over Inf).

         IF( IINFO.GT.N_SUB .AND. INFO.EQ.0 ) THEN
            INFO = 2*IOFFSET + IINFO
         ELSE IF( IINFO.LE.N_SUB .AND. IINFO.GT.0 ) THEN
            INFO = IOFFSET + IINFO
         END IF

      ELSE

         // Compute the return values for blocked code.

         // Set the number of factorized columns if the unblocked routine
         // was not called.

            K = JMAX

         // If there exits a residual matrix after the blocked code:
            // 1) compute the values of MAXC2NRMK, RELMAXC2NRMK of the
               // residual matrix, otherwise set them to ZERO;
            // 2) Set TAU(K+1:MINMN) to ZERO.

         IF( K.LT.MINMN ) THEN
            JMAXC2NRM = K + ISAMAX( N-K, WORK( K+1 ), 1 )
            MAXC2NRMK = WORK( JMAXC2NRM )
            IF( K.EQ.0 ) THEN
               RELMAXC2NRMK = ONE
            ELSE
               RELMAXC2NRMK = MAXC2NRMK / MAXC2NRM
            END IF

            DO J = K + 1, MINMN
               TAU( J ) = ZERO
            END DO

         END IF

      // END IF( J.LE.JMAX ) THEN

      END IF

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

      RETURN

      // End of SGEQP3RK

      }
