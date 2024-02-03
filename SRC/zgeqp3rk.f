      SUBROUTINE ZGEQP3RK( M, N, NRHS, KMAX, ABSTOL, RELTOL, A, LDA, K, MAXC2NRMK, RELMAXC2NRMK, JPIV, TAU, WORK, LWORK, RWORK, IWORK, INFO )
      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, KF, KMAX, LDA, LWORK, M, N, NRHS;
      double             ABSTOL,  MAXC2NRMK, RELMAXC2NRMK, RELTOL;
      // ..
      // .. Array Arguments ..
      int                IWORK( * ), JPIV( * );
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                INB, INBMIN, IXOVER;
      const              INB = 1, INBMIN = 2, IXOVER = 3 ;
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 ;
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, DONE;
      int                IINFO, IOFFSET, IWS, J, JB, JBF, JMAXB, JMAX, JMAXC2NRM, KP1, LWKOPT, MINMN, N_SUB, NB, NBMIN, NX;
      double             EPS, HUGEVAL, MAXC2NRM, SAFMIN;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLAQP2RK, ZLAQP3RK, XERBLA
      // ..
      // .. External Functions ..
      bool               DISNAN;
      int                IDAMAX, ILAENV;
      double             DLAMCH, DZNRM2;
      // EXTERNAL DISNAN, DLAMCH, DZNRM2, IDAMAX, ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test input arguments
      // ====================

      INFO = 0
      LQUERY = ( LWORK == -1 )
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( KMAX.LT.0 ) {
         INFO = -4
      } else if ( DISNAN( ABSTOL ) ) {
         INFO = -5
      } else if ( DISNAN( RELTOL ) ) {
         INFO = -6
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -8
      }

      // If the input parameters M, N, NRHS, KMAX, LDA are valid:
        // a) Test the input workspace size LWORK for the minimum
           // size requirement IWS.
        // b) Determine the optimal block size NB and optimal
           // workspace size LWKOPT to be returned in WORK(1)
           // in case of (1) LWORK < IWS, (2) LQUERY = true ,
           // (3) when routine exits.
      // Here, IWS is the miminum workspace required for unblocked
      // code.

      if ( INFO == 0 ) {
         MINMN = MIN( M, N )
         if ( MINMN == 0 ) {
            IWS = 1
            LWKOPT = 1
         } else {

            // Minimal workspace size in case of using only unblocked
            // BLAS 2 code in ZLAQP2RK.
            // 1) ZLAQP2RK: N+NRHS-1 to use in WORK array that is used
               // in ZLARF subroutine inside ZLAQP2RK to apply an
               // elementary reflector from the left.
            // TOTAL_WORK_SIZE = 3*N + NRHS - 1

            IWS = N + NRHS - 1

            // Assign to NB optimal block size.

            NB = ILAENV( INB, 'ZGEQP3RK', ' ', M, N, -1, -1 )

            // A formula for the optimal workspace size in case of using
            // both unblocked BLAS 2 in ZLAQP2RK and blocked BLAS 3 code
            // in ZLAQP3RK.
            // 1) ZGEQP3RK, ZLAQP2RK, ZLAQP3RK: 2*N to store full and
               // partial column 2-norms.
            // 2) ZLAQP2RK: N+NRHS-1 to use in WORK array that is used
               // in ZLARF subroutine to apply an elementary reflector
               // from the left.
            // 3) ZLAQP3RK: NB*(N+NRHS) to use in the work array F that
               // is used to apply a block reflector from
               // the left.
            // 4) ZLAQP3RK: NB to use in the auxilixary array AUX.
            // Sizes (2) and ((3) + (4)) should intersect, therefore
            // TOTAL_WORK_SIZE = 2*N + NB*( N+NRHS+1 ), given NBMIN=2.

            LWKOPT = 2*N + NB*( N+NRHS+1 )
         }
         WORK( 1 ) = DCMPLX( LWKOPT )

         if ( ( LWORK.LT.IWS ) && .NOT.LQUERY ) {
            INFO = -15
         }
      }

       // NOTE: The optimal workspace size is returned in WORK(1), if
             // the input parameters M, N, NRHS, KMAX, LDA are valid.

      if ( INFO != 0 ) {
         xerbla('ZGEQP3RK', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible for M=0 or N=0.

      if ( MINMN == 0 ) {
         K = 0
         MAXC2NRMK = ZERO
         RELMAXC2NRMK = ZERO
         WORK( 1 ) = DCMPLX( LWKOPT )
         RETURN
      }

      // ==================================================================

      // Initialize column pivot array JPIV.

      for (J = 1; J <= N; J++) {
         JPIV( J ) = J
      }

      // ==================================================================

      // Initialize storage for partial and exact column 2-norms.
      // a) The elements WORK(1:N) are used to store partial column
         // 2-norms of the matrix A, and may decrease in each computation
         // step; initialize to the values of complete columns 2-norms.
      // b) The elements WORK(N+1:2*N) are used to store complete column
         // 2-norms of the matrix A, they are not changed during the
         // computation; initialize the values of complete columns 2-norms.

      for (J = 1; J <= N; J++) {
         RWORK( J ) = DZNRM2( M, A( 1, J ), 1 )
         RWORK( N+J ) = RWORK( J )
      }

      // ==================================================================

      // Compute the pivot column index and the maximum column 2-norm
      // for the whole original matrix stored in A(1:M,1:N).

      KP1 = IDAMAX( N, RWORK( 1 ), 1 )

      // ==================================================================.

      if ( DISNAN( MAXC2NRM ) ) {

         // Check if the matrix A contains NaN, set INFO parameter
         // to the column number where the first NaN is found and return
         // from the routine.

         K = 0
         INFO = KP1

         // Set MAXC2NRMK and  RELMAXC2NRMK to NaN.

         MAXC2NRMK = MAXC2NRM
         RELMAXC2NRMK = MAXC2NRM

         // Array TAU is not set and contains undefined elements.

         WORK( 1 ) = DCMPLX( LWKOPT )
         RETURN
      }

      // ===================================================================

      if ( MAXC2NRM == ZERO ) {

         // Check is the matrix A is a zero matrix, set array TAU and
         // return from the routine.

         K = 0
         MAXC2NRMK = ZERO
         RELMAXC2NRMK = ZERO

         for (J = 1; J <= MINMN; J++) {
            TAU( J ) = CZERO
         }

         WORK( 1 ) = DCMPLX( LWKOPT )
         RETURN

      }

      // ===================================================================

      HUGEVAL = DLAMCH( 'Overflow' )

      if ( MAXC2NRM.GT.HUGEVAL ) {

         // Check if the matrix A contains +Inf or -Inf, set INFO parameter
         // to the column number, where the first +/-Inf  is found plus N,
         // and continue the computation.

         INFO = N + KP1

      }

      // ==================================================================

      // Quick return if possible for the case when the first
      // stopping criterion is satisfied, i.e. KMAX = 0.

      if ( KMAX == 0 ) {
         K = 0
         MAXC2NRMK = MAXC2NRM
         RELMAXC2NRMK = ONE
         for (J = 1; J <= MINMN; J++) {
            TAU( J ) = CZERO
         }
         WORK( 1 ) = DCMPLX( LWKOPT )
         RETURN
      }

      // ==================================================================

      EPS = DLAMCH('Epsilon')

      // Adjust ABSTOL

      if ( ABSTOL.GE.ZERO ) {
         SAFMIN = DLAMCH('Safe minimum')
         ABSTOL = MAX( ABSTOL, TWO*SAFMIN )
      }

      // Adjust RELTOL

      if ( RELTOL.GE.ZERO ) {
         RELTOL = MAX( RELTOL, EPS )
      }

      // ===================================================================

      // JMAX is the maximum index of the column to be factorized,
      // which is also limited by the first stopping criterion KMAX.

      JMAX = MIN( KMAX, MINMN )

      // ===================================================================

      // Quick return if possible for the case when the second or third
      // stopping criterion for the whole original matrix is satified,
      // i.e. MAXC2NRM <= ABSTOL or RELMAXC2NRM <= RELTOL
      // (which is ONE <= RELTOL).

      if ( MAXC2NRM.LE.ABSTOL .OR. ONE.LE.RELTOL ) {

         K = 0
         MAXC2NRMK = MAXC2NRM
         RELMAXC2NRMK = ONE

         for (J = 1; J <= MINMN; J++) {
            TAU( J ) = CZERO
         }

         WORK( 1 ) = DCMPLX( LWKOPT )
         RETURN
      }

      // ==================================================================
      // Factorize columns
      // ==================================================================

      // Determine the block size.

      NBMIN = 2
      NX = 0

      if ( ( NB.GT.1 ) && ( NB.LT.MINMN ) ) {

         // Determine when to cross over from blocked to unblocked code.
         // (for N less than NX, unblocked code should be used).

         NX = MAX( 0, ILAENV( IXOVER, 'ZGEQP3RK', ' ', M, N, -1, -1 ) )

         if ( NX.LT.MINMN ) {

            // Determine if workspace is large enough for blocked code.

            if ( LWORK.LT.LWKOPT ) {

               // Not enough workspace to use optimal block size that
               // is currently stored in NB.
               // Reduce NB and determine the minimum value of NB.

               NB = ( LWORK-2*N ) / ( N+1 )
               NBMIN = MAX( 2, ILAENV( INBMIN, 'ZGEQP3RK', ' ', M, N, -1, -1 ) )

            }
         }
      }

      // ==================================================================

      // DONE is the boolean flag to rerpresent the case when the
      // factorization completed in the block factorization routine,
      // before the end of the block.

      DONE = false;

      // J is the column index.

      J = 1

      // (1) Use blocked code initially.

      // JMAXB is the maximum column index of the block, when the
      // blocked code is used, is also limited by the first stopping
      // criterion KMAX.

      JMAXB = MIN( KMAX, MINMN - NX )

      if ( NB.GE.NBMIN && NB.LT.JMAX && JMAXB.GT.0 ) {

         // Loop over the column blocks of the matrix A(1:M,1:JMAXB). Here:
         // J   is the column index of a column block;
         // JB  is the column block size to pass to block factorization
             // routine in a loop step;
         // JBF is the number of columns that were actually factorized
             // that was returned by the block factorization routine
             // in a loop step, JBF <= JB;
         // N_SUB is the number of columns in the submatrix;
         // IOFFSET is the number of rows that should not be factorized.

         DO WHILE( J.LE.JMAXB )

            JB = MIN( NB, JMAXB-J+1 )
            N_SUB = N-J+1
            IOFFSET = J-1

            // Factorize JB columns among the columns A(J:N).

            zlaqp3rk(M, N_SUB, NRHS, IOFFSET, JB, ABSTOL, RELTOL, KP1, MAXC2NRM, A( 1, J ), LDA, DONE, JBF, MAXC2NRMK, RELMAXC2NRMK, JPIV( J ), TAU( J ), RWORK( J ), RWORK( N+J ), WORK( 1 ), WORK( JB+1 ), N+NRHS-J+1, IWORK, IINFO );

            // Set INFO on the first occurence of Inf.

            if ( IINFO.GT.N_SUB && INFO == 0 ) {
               INFO = 2*IOFFSET + IINFO
            }

            if ( DONE ) {

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

               if ( IINFO.LE.N_SUB && IINFO.GT.0 ) {
                  INFO = IOFFSET + IINFO
               }

               // Return from the routine.

               WORK( 1 ) = DCMPLX( LWKOPT )

               RETURN

            }

            J = J + JBF

         }

      }

      // Use unblocked code to factor the last or only block.
      // J = JMAX+1 means we factorized the maximum possible number of
      // columns, that is in ELSE clause we need to compute
      // the MAXC2NORM and RELMAXC2NORM to return after we processed
      // the blocks.

      if ( J.LE.JMAX ) {

         // N_SUB is the number of columns in the submatrix;
         // IOFFSET is the number of rows that should not be factorized.

         N_SUB = N-J+1
         IOFFSET = J-1

         zlaqp2rk(M, N_SUB, NRHS, IOFFSET, JMAX-J+1, ABSTOL, RELTOL, KP1, MAXC2NRM, A( 1, J ), LDA, KF, MAXC2NRMK, RELMAXC2NRMK, JPIV( J ), TAU( J ), RWORK( J ), RWORK( N+J ), WORK( 1 ), IINFO );

         // ABSTOL or RELTOL criterion is satisfied when the number of
         // the factorized columns KF is smaller then the  number
         // of columns JMAX-J+1 supplied to be factorized by the
         // unblocked routine, we can return from
         // the routine. Perform the following before returning:
            // a) Set the number of factorized columns K,
            // b) MAXC2NRMK and RELMAXC2NRMK are returned by the
               // unblocked factorization routine above.

         K = J - 1 + KF

         // Set INFO on the first exception occurence.

         // Set INFO on the first exception occurence of Inf or NaN,
         // (NaN takes precedence over Inf).

         if ( IINFO.GT.N_SUB && INFO == 0 ) {
            INFO = 2*IOFFSET + IINFO
         } else if ( IINFO.LE.N_SUB && IINFO.GT.0 ) {
            INFO = IOFFSET + IINFO
         }

      } else {

         // Compute the return values for blocked code.

         // Set the number of factorized columns if the unblocked routine
         // was not called.

            K = JMAX

         // If there exits a residual matrix after the blocked code:
            // 1) compute the values of MAXC2NRMK, RELMAXC2NRMK of the
               // residual matrix, otherwise set them to ZERO;
            // 2) Set TAU(K+1:MINMN) to ZERO.

         if ( K.LT.MINMN ) {
            JMAXC2NRM = K + IDAMAX( N-K, RWORK( K+1 ), 1 )
            MAXC2NRMK = RWORK( JMAXC2NRM )
            if ( K == 0 ) {
               RELMAXC2NRMK = ONE
            } else {
               RELMAXC2NRMK = MAXC2NRMK / MAXC2NRM
            }

            for (J = K + 1; J <= MINMN; J++) {
               TAU( J ) = CZERO
            }

         } else {
            MAXC2NRMK = ZERO
            RELMAXC2NRMK = ZERO

         }

      // END IF( J.LE.JMAX ) THEN

      }

      WORK( 1 ) = DCMPLX( LWKOPT )

      RETURN

      // End of ZGEQP3RK

      }
