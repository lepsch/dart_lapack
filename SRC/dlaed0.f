      SUBROUTINE DLAED0( ICOMPQ, QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ICOMPQ, INFO, LDQ, LDQS, N, QSIZ;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             D( * ), E( * ), Q( LDQ, * ), QSTORE( LDQS, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      PARAMETER          ( ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0 )
      // ..
      // .. Local Scalars ..
      int                CURLVL, CURPRB, CURR, I, IGIVCL, IGIVNM, IGIVPT, INDXQ, IPERM, IPRMPT, IQ, IQPTR, IWREM, J, K, LGN, MATSIZ, MSD2, SMLSIZ, SMM1, SPM1, SPM2, SUBMAT, SUBPBS, TLVLS;
      double             TEMP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMM, DLACPY, DLAED1, DLAED7, DSTEQR, XERBLA
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

      IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.2 ) THEN
         INFO = -1
      ELSE IF( ( ICOMPQ.EQ.1 ) .AND. ( QSIZ.LT.MAX( 0, N ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDQS.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED0', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      SMLSIZ = ILAENV( 9, 'DLAED0', ' ', 0, 0, 0, 0 )

      // Determine the size and placement of the submatrices, and save in
     t // he leading elements of IWORK.

      IWORK( 1 ) = N
      SUBPBS = 1
      TLVLS = 0
   10 CONTINUE
      IF( IWORK( SUBPBS ).GT.SMLSIZ ) THEN
         DO 20 J = SUBPBS, 1, -1
            IWORK( 2*J ) = ( IWORK( J )+1 ) / 2
            IWORK( 2*J-1 ) = IWORK( J ) / 2
   20    CONTINUE
         TLVLS = TLVLS + 1
         SUBPBS = 2*SUBPBS
         GO TO 10
      END IF
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
      IF( ICOMPQ.NE.2 ) THEN

         // Set up workspaces for eigenvalues only/accumulate new vectors
         // routine

         TEMP = LOG( DBLE( N ) ) / LOG( TWO )
         LGN = INT( TEMP )
         IF( 2**LGN.LT.N ) LGN = LGN + 1          IF( 2**LGN.LT.N ) LGN = LGN + 1
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
   50    CONTINUE
         IWORK( IQPTR ) = 1
      END IF

      // Solve each submatrix eigenproblem at the bottom of the divide and
      // conquer tree.

      CURR = 0
      DO 70 I = 0, SPM1
         IF( I.EQ.0 ) THEN
            SUBMAT = 1
            MATSIZ = IWORK( 1 )
         ELSE
            SUBMAT = IWORK( I ) + 1
            MATSIZ = IWORK( I+1 ) - IWORK( I )
         END IF
         IF( ICOMPQ.EQ.2 ) THEN
            CALL DSTEQR( 'I', MATSIZ, D( SUBMAT ), E( SUBMAT ), Q( SUBMAT, SUBMAT ), LDQ, WORK, INFO )             IF( INFO.NE.0 ) GO TO 130
         ELSE
            CALL DSTEQR( 'I', MATSIZ, D( SUBMAT ), E( SUBMAT ), WORK( IQ-1+IWORK( IQPTR+CURR ) ), MATSIZ, WORK, INFO )
            IF( INFO.NE.0 ) GO TO 130
            IF( ICOMPQ.EQ.1 ) THEN
               CALL DGEMM( 'N', 'N', QSIZ, MATSIZ, MATSIZ, ONE, Q( 1, SUBMAT ), LDQ, WORK( IQ-1+IWORK( IQPTR+ CURR ) ), MATSIZ, ZERO, QSTORE( 1, SUBMAT ), LDQS )
            END IF
            IWORK( IQPTR+CURR+1 ) = IWORK( IQPTR+CURR ) + MATSIZ**2
            CURR = CURR + 1
         END IF
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
      IF( SUBPBS.GT.1 ) THEN
         SPM2 = SUBPBS - 2
         DO 90 I = 0, SPM2, 2
            IF( I.EQ.0 ) THEN
               SUBMAT = 1
               MATSIZ = IWORK( 2 )
               MSD2 = IWORK( 1 )
               CURPRB = 0
            ELSE
               SUBMAT = IWORK( I ) + 1
               MATSIZ = IWORK( I+2 ) - IWORK( I )
               MSD2 = MATSIZ / 2
               CURPRB = CURPRB + 1
            END IF

      // Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
      // into an eigensystem of size MATSIZ.
      // DLAED1 is used only for the full eigensystem of a tridiagonal
      // matrix.
      // DLAED7 handles the cases in which eigenvalues only or eigenvalues
      // and eigenvectors of a full symmetric matrix (which was reduced to
     t // ridiagonal form) are desired.

            IF( ICOMPQ.EQ.2 ) THEN
               CALL DLAED1( MATSIZ, D( SUBMAT ), Q( SUBMAT, SUBMAT ), LDQ, IWORK( INDXQ+SUBMAT ), E( SUBMAT+MSD2-1 ), MSD2, WORK, IWORK( SUBPBS+1 ), INFO )
            ELSE
               CALL DLAED7( ICOMPQ, MATSIZ, QSIZ, TLVLS, CURLVL, CURPRB, D( SUBMAT ), QSTORE( 1, SUBMAT ), LDQS, IWORK( INDXQ+SUBMAT ), E( SUBMAT+MSD2-1 ), MSD2, WORK( IQ ), IWORK( IQPTR ), IWORK( IPRMPT ), IWORK( IPERM ), IWORK( IGIVPT ), IWORK( IGIVCL ), WORK( IGIVNM ), WORK( IWREM ), IWORK( SUBPBS+1 ), INFO )
            END IF
            IF( INFO.NE.0 ) GO TO 130
            IWORK( I / 2+1 ) = IWORK( I+2 )
   90    CONTINUE
         SUBPBS = SUBPBS / 2
         CURLVL = CURLVL + 1
         GO TO 80
      END IF

      // end while

      // Re-merge the eigenvalues/vectors which were deflated at the final
      // merge step.

      IF( ICOMPQ.EQ.1 ) THEN
         DO 100 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
            CALL DCOPY( QSIZ, QSTORE( 1, J ), 1, Q( 1, I ), 1 )
  100    CONTINUE
         CALL DCOPY( N, WORK, 1, D, 1 )
      ELSE IF( ICOMPQ.EQ.2 ) THEN
         DO 110 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
            CALL DCOPY( N, Q( 1, J ), 1, WORK( N*I+1 ), 1 )
  110    CONTINUE
         CALL DCOPY( N, WORK, 1, D, 1 )
         CALL DLACPY( 'A', N, N, WORK( N+1 ), N, Q, LDQ )
      ELSE
         DO 120 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
  120    CONTINUE
         CALL DCOPY( N, WORK, 1, D, 1 )
      END IF
      GO TO 140

  130 CONTINUE
      INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1

  140 CONTINUE
      RETURN

      // End of DLAED0

      END
