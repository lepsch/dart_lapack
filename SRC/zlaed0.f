      SUBROUTINE ZLAED0( QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, RWORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDQ, LDQS, N, QSIZ;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             D( * ), E( * ), RWORK( * );
      COMPLEX*16         Q( LDQ, * ), QSTORE( LDQS, * )
      // ..

*  =====================================================================

*  Warning:      N could be as big as QSIZ!

      // .. Parameters ..
      double             TWO;
      const              TWO = 2.D+0 ;
      // ..
      // .. Local Scalars ..
      int                CURLVL, CURPRB, CURR, I, IGIVCL, IGIVNM, IGIVPT, INDXQ, IPERM, IPRMPT, IQ, IQPTR, IWREM, J, K, LGN, LL, MATSIZ, MSD2, SMLSIZ, SMM1, SPM1, SPM2, SUBMAT, SUBPBS, TLVLS;
      double             TEMP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DSTEQR, XERBLA, ZCOPY, ZLACRM, ZLAED7
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, INT, LOG, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      // IF( ICOMPQ .LT. 0 .OR. ICOMPQ .GT. 2 ) THEN
         // INFO = -1
      // ELSE IF( ( ICOMPQ .EQ. 1 ) .AND. ( QSIZ .LT. MAX( 0, N ) ) )
*    $        THEN
      if ( QSIZ.LT.MAX( 0, N ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDQ.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDQS.LT.MAX( 1, N ) ) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         xerbla('ZLAED0', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      SMLSIZ = ILAENV( 9, 'ZLAED0', ' ', 0, 0, 0, 0 )

      // Determine the size and placement of the submatrices, and save in
      // the leading elements of IWORK.

      IWORK( 1 ) = N
      SUBPBS = 1
      TLVLS = 0
   10 CONTINUE
      if ( IWORK( SUBPBS ).GT.SMLSIZ ) {
         DO 20 J = SUBPBS, 1, -1
            IWORK( 2*J ) = ( IWORK( J )+1 ) / 2
            IWORK( 2*J-1 ) = IWORK( J ) / 2
   20    CONTINUE
         TLVLS = TLVLS + 1
         SUBPBS = 2*SUBPBS
         GO TO 10
      }
      DO 30 J = 2, SUBPBS
         IWORK( J ) = IWORK( J ) + IWORK( J-1 )
   30 CONTINUE

      // Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
      // using rank-1 modifications (cuts).

      SPM1 = SUBPBS - 1
      DO 40 I = 1, SPM1
         SUBMAT = IWORK( I ) + 1
         SMM1 = SUBMAT - 1
         D( SMM1 ) = D( SMM1 ) - ABS( E( SMM1 ) )
         D( SUBMAT ) = D( SUBMAT ) - ABS( E( SMM1 ) )
   40 CONTINUE

      INDXQ = 4*N + 3

      // Set up workspaces for eigenvalues only/accumulate new vectors
      // routine

      TEMP = LOG( DBLE( N ) ) / LOG( TWO )
      LGN = INT( TEMP )
      IF( 2**LGN.LT.N ) LGN = LGN + 1       IF( 2**LGN.LT.N ) LGN = LGN + 1
      IPRMPT = INDXQ + N + 1
      IPERM = IPRMPT + N*LGN
      IQPTR = IPERM + N*LGN
      IGIVPT = IQPTR + N + 2
      IGIVCL = IGIVPT + N*LGN

      IGIVNM = 1
      IQ = IGIVNM + 2*N*LGN
      IWREM = IQ + N**2 + 1
      // Initialize pointers
      DO 50 I = 0, SUBPBS
         IWORK( IPRMPT+I ) = 1
         IWORK( IGIVPT+I ) = 1
   50 CONTINUE
      IWORK( IQPTR ) = 1

      // Solve each submatrix eigenproblem at the bottom of the divide and
      // conquer tree.

      CURR = 0
      DO 70 I = 0, SPM1
         if ( I.EQ.0 ) {
            SUBMAT = 1
            MATSIZ = IWORK( 1 )
         } else {
            SUBMAT = IWORK( I ) + 1
            MATSIZ = IWORK( I+1 ) - IWORK( I )
         }
         LL = IQ - 1 + IWORK( IQPTR+CURR )
         dsteqr('I', MATSIZ, D( SUBMAT ), E( SUBMAT ), RWORK( LL ), MATSIZ, RWORK, INFO )          CALL ZLACRM( QSIZ, MATSIZ, Q( 1, SUBMAT ), LDQ, RWORK( LL ), MATSIZ, QSTORE( 1, SUBMAT ), LDQS, RWORK( IWREM ) );
         IWORK( IQPTR+CURR+1 ) = IWORK( IQPTR+CURR ) + MATSIZ**2
         CURR = CURR + 1
         if ( INFO.GT.0 ) {
            INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1
            RETURN
         }
         K = 1
         DO 60 J = SUBMAT, IWORK( I+1 )
            IWORK( INDXQ+J ) = K
            K = K + 1
   60    CONTINUE
   70 CONTINUE

      // Successively merge eigensystems of adjacent submatrices
      // into eigensystem for the corresponding larger matrix.

      // while ( SUBPBS > 1 )

      CURLVL = 1
   80 CONTINUE
      if ( SUBPBS.GT.1 ) {
         SPM2 = SUBPBS - 2
         DO 90 I = 0, SPM2, 2
            if ( I.EQ.0 ) {
               SUBMAT = 1
               MATSIZ = IWORK( 2 )
               MSD2 = IWORK( 1 )
               CURPRB = 0
            } else {
               SUBMAT = IWORK( I ) + 1
               MATSIZ = IWORK( I+2 ) - IWORK( I )
               MSD2 = MATSIZ / 2
               CURPRB = CURPRB + 1
            }

      // Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
      // into an eigensystem of size MATSIZ.  ZLAED7 handles the case
      // when the eigenvectors of a full or band Hermitian matrix (which
      // was reduced to tridiagonal form) are desired.

      // I am free to use Q as a valuable working space until Loop 150.

            zlaed7(MATSIZ, MSD2, QSIZ, TLVLS, CURLVL, CURPRB, D( SUBMAT ), QSTORE( 1, SUBMAT ), LDQS, E( SUBMAT+MSD2-1 ), IWORK( INDXQ+SUBMAT ), RWORK( IQ ), IWORK( IQPTR ), IWORK( IPRMPT ), IWORK( IPERM ), IWORK( IGIVPT ), IWORK( IGIVCL ), RWORK( IGIVNM ), Q( 1, SUBMAT ), RWORK( IWREM ), IWORK( SUBPBS+1 ), INFO );
            if ( INFO.GT.0 ) {
               INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1
               RETURN
            }
            IWORK( I / 2+1 ) = IWORK( I+2 )
   90    CONTINUE
         SUBPBS = SUBPBS / 2
         CURLVL = CURLVL + 1
         GO TO 80
      }

      // end while

      // Re-merge the eigenvalues/vectors which were deflated at the final
      // merge step.

      DO 100 I = 1, N
         J = IWORK( INDXQ+I )
         RWORK( I ) = D( J )
         zcopy(QSIZ, QSTORE( 1, J ), 1, Q( 1, I ), 1 );
  100 CONTINUE
      dcopy(N, RWORK, 1, D, 1 );

      RETURN

      // End of ZLAED0

      }
