      SUBROUTINE SDRGSX( NSIZE, NCMAX, THRESH, NIN, NOUT, A, LDA, B, AI, BI, Z, Q, ALPHAR, ALPHAI, BETA, C, LDC, S, WORK, LWORK, IWORK, LIWORK, BWORK, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, LDC, LIWORK, LWORK, NCMAX, NIN, NOUT, NSIZE
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            BWORK( * )
      int                IWORK( * )
      REAL               A( LDA, * ), AI( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDA, * ), BETA( * ), BI( LDA, * ), C( LDC, * ), Q( LDA, * ), S( * ), WORK( * ), Z( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TEN
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TEN = 1.0E+1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILABAD
      String             SENSE;
      int                BDSPAC, I, I1, IFUNC, IINFO, J, LINFO, MAXWRK, MINWRK, MM, MN2, NERRS, NPTKNT, NTEST, NTESTT, PRTYPE, QBA, QBB
      REAL               ABNRM, BIGNUM, DIFTRU, PLTRU, SMLNUM, TEMP1, TEMP2, THRSH2, ULP, ULPINV, WEIGHT
*     ..
*     .. Local Arrays ..
      REAL               DIFEST( 2 ), PL( 2 ), RESULT( 10 )
*     ..
*     .. External Functions ..
      LOGICAL            SLCTSX
      int                ILAENV
      REAL               SLAMCH, SLANGE
      EXTERNAL           SLCTSX, ILAENV, SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALASVM, SGESVD, SGET51, SGET53, SGGESX, SLACPY, SLAKF2, SLASET, SLATM5, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Scalars in Common ..
      LOGICAL            FS
      int                K, M, MPLUSN, N
*     ..
*     .. Common blocks ..
      COMMON             / MN / M, N, MPLUSN, K, FS
*     ..
*     .. Executable Statements ..
*
*     Check for errors
*
      IF( NSIZE.LT.0 ) THEN
         INFO = -1
      ELSE IF( THRESH.LT.ZERO ) THEN
         INFO = -2
      ELSE IF( NIN.LE.0 ) THEN
         INFO = -3
      ELSE IF( NOUT.LE.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.1 .OR. LDA.LT.NSIZE ) THEN
         INFO = -6
      ELSE IF( LDC.LT.1 .OR. LDC.LT.NSIZE*NSIZE / 2 ) THEN
         INFO = -17
      ELSE IF( LIWORK.LT.NSIZE+6 ) THEN
         INFO = -21
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV.)
*
      MINWRK = 1
      IF( INFO.EQ.0 .AND. LWORK.GE.1 ) THEN
c        MINWRK = MAX( 10*( NSIZE+1 ), 5*NSIZE*NSIZE / 2-2 )
         MINWRK = MAX( 10*( NSIZE+1 ), 5*NSIZE*NSIZE / 2 )
*
*        workspace for sggesx
*
         MAXWRK = 9*( NSIZE+1 ) + NSIZE* ILAENV( 1, 'SGEQRF', ' ', NSIZE, 1, NSIZE, 0 )          MAXWRK = MAX( MAXWRK, 9*( NSIZE+1 )+NSIZE* ILAENV( 1, 'SORGQR', ' ', NSIZE, 1, NSIZE, -1 ) )
*
*        workspace for sgesvd
*
         BDSPAC = 5*NSIZE*NSIZE / 2
         MAXWRK = MAX( MAXWRK, 3*NSIZE*NSIZE / 2+NSIZE*NSIZE* ILAENV( 1, 'SGEBRD', ' ', NSIZE*NSIZE / 2, NSIZE*NSIZE / 2, -1, -1 ) )
         MAXWRK = MAX( MAXWRK, BDSPAC )
*
         MAXWRK = MAX( MAXWRK, MINWRK )
*
         WORK( 1 ) = MAXWRK
      END IF
*
      IF( LWORK.LT.MINWRK ) INFO = -19
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SDRGSX', -INFO )
         RETURN
      END IF
*
*     Important constants
*
      ULP = SLAMCH( 'P' )
      ULPINV = ONE / ULP
      SMLNUM = SLAMCH( 'S' ) / ULP
      BIGNUM = ONE / SMLNUM
      THRSH2 = TEN*THRESH
      NTESTT = 0
      NERRS = 0
*
*     Go to the tests for read-in matrix pairs
*
      IFUNC = 0
      IF( NSIZE.EQ.0 ) GO TO 70
*
*     Test the built-in matrix pairs.
*     Loop over different functions (IFUNC) of SGGESX, types (PRTYPE)
*     of test matrices, different size (M+N)
*
      PRTYPE = 0
      QBA = 3
      QBB = 4
      WEIGHT = SQRT( ULP )
*
      DO 60 IFUNC = 0, 3
         DO 50 PRTYPE = 1, 5
            DO 40 M = 1, NSIZE - 1
               DO 30 N = 1, NSIZE - M
*
                  WEIGHT = ONE / WEIGHT
                  MPLUSN = M + N
*
*                 Generate test matrices
*
                  FS = .TRUE.
                  K = 0
*
                  CALL SLASET( 'Full', MPLUSN, MPLUSN, ZERO, ZERO, AI, LDA )                   CALL SLASET( 'Full', MPLUSN, MPLUSN, ZERO, ZERO, BI, LDA )
*
                  CALL SLATM5( PRTYPE, M, N, AI, LDA, AI( M+1, M+1 ), LDA, AI( 1, M+1 ), LDA, BI, LDA, BI( M+1, M+1 ), LDA, BI( 1, M+1 ), LDA, Q, LDA, Z, LDA, WEIGHT, QBA, QBB )
*
*                 Compute the Schur factorization and swapping the
*                 m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.
*                 Swapping is accomplished via the function SLCTSX
*                 which is supplied below.
*
                  IF( IFUNC.EQ.0 ) THEN
                     SENSE = 'N'
                  ELSE IF( IFUNC.EQ.1 ) THEN
                     SENSE = 'E'
                  ELSE IF( IFUNC.EQ.2 ) THEN
                     SENSE = 'V'
                  ELSE IF( IFUNC.EQ.3 ) THEN
                     SENSE = 'B'
                  END IF
*
                  CALL SLACPY( 'Full', MPLUSN, MPLUSN, AI, LDA, A, LDA )
                  CALL SLACPY( 'Full', MPLUSN, MPLUSN, BI, LDA, B, LDA )
*
                  CALL SGGESX( 'V', 'V', 'S', SLCTSX, SENSE, MPLUSN, AI, LDA, BI, LDA, MM, ALPHAR, ALPHAI, BETA, Q, LDA, Z, LDA, PL, DIFEST, WORK, LWORK, IWORK, LIWORK, BWORK, LINFO )
*
                  IF( LINFO.NE.0 .AND. LINFO.NE.MPLUSN+2 ) THEN
                     RESULT( 1 ) = ULPINV
                     WRITE( NOUT, FMT = 9999 )'SGGESX', LINFO, MPLUSN, PRTYPE
                     INFO = LINFO
                     GO TO 30
                  END IF
*
*                 Compute the norm(A, B)
*
                  CALL SLACPY( 'Full', MPLUSN, MPLUSN, AI, LDA, WORK, MPLUSN )                   CALL SLACPY( 'Full', MPLUSN, MPLUSN, BI, LDA, WORK( MPLUSN*MPLUSN+1 ), MPLUSN )                   ABNRM = SLANGE( 'Fro', MPLUSN, 2*MPLUSN, WORK, MPLUSN, WORK )
*
*                 Do tests (1) to (4)
*
                  CALL SGET51( 1, MPLUSN, A, LDA, AI, LDA, Q, LDA, Z, LDA, WORK, RESULT( 1 ) )                   CALL SGET51( 1, MPLUSN, B, LDA, BI, LDA, Q, LDA, Z, LDA, WORK, RESULT( 2 ) )                   CALL SGET51( 3, MPLUSN, B, LDA, BI, LDA, Q, LDA, Q, LDA, WORK, RESULT( 3 ) )                   CALL SGET51( 3, MPLUSN, B, LDA, BI, LDA, Z, LDA, Z, LDA, WORK, RESULT( 4 ) )
                  NTEST = 4
*
*                 Do tests (5) and (6): check Schur form of A and
*                 compare eigenvalues with diagonals.
*
                  TEMP1 = ZERO
                  RESULT( 5 ) = ZERO
                  RESULT( 6 ) = ZERO
*
                  DO 10 J = 1, MPLUSN
                     ILABAD = .FALSE.
                     IF( ALPHAI( J ).EQ.ZERO ) THEN
                        TEMP2 = ( ABS( ALPHAR( J )-AI( J, J ) ) / MAX( SMLNUM, ABS( ALPHAR( J ) ), ABS( AI( J, J ) ) )+ ABS( BETA( J )-BI( J, J ) ) / MAX( SMLNUM, ABS( BETA( J ) ), ABS( BI( J, J ) ) ) ) / ULP
                        IF( J.LT.MPLUSN ) THEN
                           IF( AI( J+1, J ).NE.ZERO ) THEN
                              ILABAD = .TRUE.
                              RESULT( 5 ) = ULPINV
                           END IF
                        END IF
                        IF( J.GT.1 ) THEN
                           IF( AI( J, J-1 ).NE.ZERO ) THEN
                              ILABAD = .TRUE.
                              RESULT( 5 ) = ULPINV
                           END IF
                        END IF
                     ELSE
                        IF( ALPHAI( J ).GT.ZERO ) THEN
                           I1 = J
                        ELSE
                           I1 = J - 1
                        END IF
                        IF( I1.LE.0 .OR. I1.GE.MPLUSN ) THEN
                           ILABAD = .TRUE.
                        ELSE IF( I1.LT.MPLUSN-1 ) THEN
                           IF( AI( I1+2, I1+1 ).NE.ZERO ) THEN
                              ILABAD = .TRUE.
                              RESULT( 5 ) = ULPINV
                           END IF
                        ELSE IF( I1.GT.1 ) THEN
                           IF( AI( I1, I1-1 ).NE.ZERO ) THEN
                              ILABAD = .TRUE.
                              RESULT( 5 ) = ULPINV
                           END IF
                        END IF
                        IF( .NOT.ILABAD ) THEN
                           CALL SGET53( AI( I1, I1 ), LDA, BI( I1, I1 ), LDA, BETA( J ), ALPHAR( J ), ALPHAI( J ), TEMP2, IINFO )
                           IF( IINFO.GE.3 ) THEN
                              WRITE( NOUT, FMT = 9997 )IINFO, J, MPLUSN, PRTYPE
                              INFO = ABS( IINFO )
                           END IF
                        ELSE
                           TEMP2 = ULPINV
                        END IF
                     END IF
                     TEMP1 = MAX( TEMP1, TEMP2 )
                     IF( ILABAD ) THEN
                        WRITE( NOUT, FMT = 9996 )J, MPLUSN, PRTYPE
                     END IF
   10             CONTINUE
                  RESULT( 6 ) = TEMP1
                  NTEST = NTEST + 2
*
*                 Test (7) (if sorting worked)
*
                  RESULT( 7 ) = ZERO
                  IF( LINFO.EQ.MPLUSN+3 ) THEN
                     RESULT( 7 ) = ULPINV
                  ELSE IF( MM.NE.N ) THEN
                     RESULT( 7 ) = ULPINV
                  END IF
                  NTEST = NTEST + 1
*
*                 Test (8): compare the estimated value DIF and its
*                 value. first, compute the exact DIF.
*
                  RESULT( 8 ) = ZERO
                  MN2 = MM*( MPLUSN-MM )*2
                  IF( IFUNC.GE.2 .AND. MN2.LE.NCMAX*NCMAX ) THEN
*
*                    Note: for either following two causes, there are
*                    almost same number of test cases fail the test.
*
                     CALL SLAKF2( MM, MPLUSN-MM, AI, LDA, AI( MM+1, MM+1 ), BI, BI( MM+1, MM+1 ), C, LDC )
*
                     CALL SGESVD( 'N', 'N', MN2, MN2, C, LDC, S, WORK, 1, WORK( 2 ), 1, WORK( 3 ), LWORK-2, INFO )
                     DIFTRU = S( MN2 )
*
                     IF( DIFEST( 2 ).EQ.ZERO ) THEN
                        IF( DIFTRU.GT.ABNRM*ULP ) RESULT( 8 ) = ULPINV
                     ELSE IF( DIFTRU.EQ.ZERO ) THEN
                        IF( DIFEST( 2 ).GT.ABNRM*ULP ) RESULT( 8 ) = ULPINV                      ELSE IF( ( DIFTRU.GT.THRSH2*DIFEST( 2 ) ) .OR. ( DIFTRU*THRSH2.LT.DIFEST( 2 ) ) ) THEN                         RESULT( 8 ) = MAX( DIFTRU / DIFEST( 2 ), DIFEST( 2 ) / DIFTRU )
                     END IF
                     NTEST = NTEST + 1
                  END IF
*
*                 Test (9)
*
                  RESULT( 9 ) = ZERO
                  IF( LINFO.EQ.( MPLUSN+2 ) ) THEN
                     IF( DIFTRU.GT.ABNRM*ULP ) RESULT( 9 ) = ULPINV                      IF( ( IFUNC.GT.1 ) .AND. ( DIFEST( 2 ).NE.ZERO ) ) RESULT( 9 ) = ULPINV                      IF( ( IFUNC.EQ.1 ) .AND. ( PL( 1 ).NE.ZERO ) ) RESULT( 9 ) = ULPINV
                     NTEST = NTEST + 1
                  END IF
*
                  NTESTT = NTESTT + NTEST
*
*                 Print out tests which fail.
*
                  DO 20 J = 1, 9
                     IF( RESULT( J ).GE.THRESH ) THEN
*
*                       If this is the first test to fail,
*                       print a header to the data file.
*
                        IF( NERRS.EQ.0 ) THEN
                           WRITE( NOUT, FMT = 9995 )'SGX'
*
*                          Matrix types
*
                           WRITE( NOUT, FMT = 9993 )
*
*                          Tests performed
*
                           WRITE( NOUT, FMT = 9992 )'orthogonal', '''', 'transpose', ( '''', I = 1, 4 )
*
                        END IF
                        NERRS = NERRS + 1
                        IF( RESULT( J ).LT.10000.0 ) THEN
                           WRITE( NOUT, FMT = 9991 )MPLUSN, PRTYPE, WEIGHT, M, J, RESULT( J )
                        ELSE
                           WRITE( NOUT, FMT = 9990 )MPLUSN, PRTYPE, WEIGHT, M, J, RESULT( J )
                        END IF
                     END IF
   20             CONTINUE
*
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
*
      GO TO 150
*
   70 CONTINUE
*
*     Read in data from file to check accuracy of condition estimation
*     Read input data until N=0
*
      NPTKNT = 0
*
   80 CONTINUE
      READ( NIN, FMT = *, END = 140 )MPLUSN
      IF( MPLUSN.EQ.0 ) GO TO 140
      READ( NIN, FMT = *, END = 140 )N
      DO 90 I = 1, MPLUSN
         READ( NIN, FMT = * )( AI( I, J ), J = 1, MPLUSN )
   90 CONTINUE
      DO 100 I = 1, MPLUSN
         READ( NIN, FMT = * )( BI( I, J ), J = 1, MPLUSN )
  100 CONTINUE
      READ( NIN, FMT = * )PLTRU, DIFTRU
*
      NPTKNT = NPTKNT + 1
      FS = .TRUE.
      K = 0
      M = MPLUSN - N
*
      CALL SLACPY( 'Full', MPLUSN, MPLUSN, AI, LDA, A, LDA )
      CALL SLACPY( 'Full', MPLUSN, MPLUSN, BI, LDA, B, LDA )
*
*     Compute the Schur factorization while swapping the
*     m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.
*
      CALL SGGESX( 'V', 'V', 'S', SLCTSX, 'B', MPLUSN, AI, LDA, BI, LDA, MM, ALPHAR, ALPHAI, BETA, Q, LDA, Z, LDA, PL, DIFEST, WORK, LWORK, IWORK, LIWORK, BWORK, LINFO )
*
      IF( LINFO.NE.0 .AND. LINFO.NE.MPLUSN+2 ) THEN
         RESULT( 1 ) = ULPINV
         WRITE( NOUT, FMT = 9998 )'SGGESX', LINFO, MPLUSN, NPTKNT
         GO TO 130
      END IF
*
*     Compute the norm(A, B)
*        (should this be norm of (A,B) or (AI,BI)?)
*
      CALL SLACPY( 'Full', MPLUSN, MPLUSN, AI, LDA, WORK, MPLUSN )
      CALL SLACPY( 'Full', MPLUSN, MPLUSN, BI, LDA, WORK( MPLUSN*MPLUSN+1 ), MPLUSN )
      ABNRM = SLANGE( 'Fro', MPLUSN, 2*MPLUSN, WORK, MPLUSN, WORK )
*
*     Do tests (1) to (4)
*
      CALL SGET51( 1, MPLUSN, A, LDA, AI, LDA, Q, LDA, Z, LDA, WORK, RESULT( 1 ) )       CALL SGET51( 1, MPLUSN, B, LDA, BI, LDA, Q, LDA, Z, LDA, WORK, RESULT( 2 ) )       CALL SGET51( 3, MPLUSN, B, LDA, BI, LDA, Q, LDA, Q, LDA, WORK, RESULT( 3 ) )       CALL SGET51( 3, MPLUSN, B, LDA, BI, LDA, Z, LDA, Z, LDA, WORK, RESULT( 4 ) )
*
*     Do tests (5) and (6): check Schur form of A and compare
*     eigenvalues with diagonals.
*
      NTEST = 6
      TEMP1 = ZERO
      RESULT( 5 ) = ZERO
      RESULT( 6 ) = ZERO
*
      DO 110 J = 1, MPLUSN
         ILABAD = .FALSE.
         IF( ALPHAI( J ).EQ.ZERO ) THEN
            TEMP2 = ( ABS( ALPHAR( J )-AI( J, J ) ) / MAX( SMLNUM, ABS( ALPHAR( J ) ), ABS( AI( J, J ) ) )+ABS( BETA( J )-BI( J, J ) ) / MAX( SMLNUM, ABS( BETA( J ) ), ABS( BI( J, J ) ) ) ) / ULP
            IF( J.LT.MPLUSN ) THEN
               IF( AI( J+1, J ).NE.ZERO ) THEN
                  ILABAD = .TRUE.
                  RESULT( 5 ) = ULPINV
               END IF
            END IF
            IF( J.GT.1 ) THEN
               IF( AI( J, J-1 ).NE.ZERO ) THEN
                  ILABAD = .TRUE.
                  RESULT( 5 ) = ULPINV
               END IF
            END IF
         ELSE
            IF( ALPHAI( J ).GT.ZERO ) THEN
               I1 = J
            ELSE
               I1 = J - 1
            END IF
            IF( I1.LE.0 .OR. I1.GE.MPLUSN ) THEN
               ILABAD = .TRUE.
            ELSE IF( I1.LT.MPLUSN-1 ) THEN
               IF( AI( I1+2, I1+1 ).NE.ZERO ) THEN
                  ILABAD = .TRUE.
                  RESULT( 5 ) = ULPINV
               END IF
            ELSE IF( I1.GT.1 ) THEN
               IF( AI( I1, I1-1 ).NE.ZERO ) THEN
                  ILABAD = .TRUE.
                  RESULT( 5 ) = ULPINV
               END IF
            END IF
            IF( .NOT.ILABAD ) THEN
               CALL SGET53( AI( I1, I1 ), LDA, BI( I1, I1 ), LDA, BETA( J ), ALPHAR( J ), ALPHAI( J ), TEMP2, IINFO )
               IF( IINFO.GE.3 ) THEN
                  WRITE( NOUT, FMT = 9997 )IINFO, J, MPLUSN, NPTKNT
                  INFO = ABS( IINFO )
               END IF
            ELSE
               TEMP2 = ULPINV
            END IF
         END IF
         TEMP1 = MAX( TEMP1, TEMP2 )
         IF( ILABAD ) THEN
            WRITE( NOUT, FMT = 9996 )J, MPLUSN, NPTKNT
         END IF
  110 CONTINUE
      RESULT( 6 ) = TEMP1
*
*     Test (7) (if sorting worked)  <--------- need to be checked.
*
      NTEST = 7
      RESULT( 7 ) = ZERO
      IF( LINFO.EQ.MPLUSN+3 ) RESULT( 7 ) = ULPINV
*
*     Test (8): compare the estimated value of DIF and its true value.
*
      NTEST = 8
      RESULT( 8 ) = ZERO
      IF( DIFEST( 2 ).EQ.ZERO ) THEN
         IF( DIFTRU.GT.ABNRM*ULP ) RESULT( 8 ) = ULPINV
      ELSE IF( DIFTRU.EQ.ZERO ) THEN
         IF( DIFEST( 2 ).GT.ABNRM*ULP ) RESULT( 8 ) = ULPINV       ELSE IF( ( DIFTRU.GT.THRSH2*DIFEST( 2 ) ) .OR. ( DIFTRU*THRSH2.LT.DIFEST( 2 ) ) ) THEN
         RESULT( 8 ) = MAX( DIFTRU / DIFEST( 2 ), DIFEST( 2 ) / DIFTRU )
      END IF
*
*     Test (9)
*
      NTEST = 9
      RESULT( 9 ) = ZERO
      IF( LINFO.EQ.( MPLUSN+2 ) ) THEN
         IF( DIFTRU.GT.ABNRM*ULP ) RESULT( 9 ) = ULPINV          IF( ( IFUNC.GT.1 ) .AND. ( DIFEST( 2 ).NE.ZERO ) ) RESULT( 9 ) = ULPINV          IF( ( IFUNC.EQ.1 ) .AND. ( PL( 1 ).NE.ZERO ) ) RESULT( 9 ) = ULPINV
      END IF
*
*     Test (10): compare the estimated value of PL and it true value.
*
      NTEST = 10
      RESULT( 10 ) = ZERO
      IF( PL( 1 ).EQ.ZERO ) THEN
         IF( PLTRU.GT.ABNRM*ULP ) RESULT( 10 ) = ULPINV
      ELSE IF( PLTRU.EQ.ZERO ) THEN
         IF( PL( 1 ).GT.ABNRM*ULP ) RESULT( 10 ) = ULPINV       ELSE IF( ( PLTRU.GT.THRESH*PL( 1 ) ) .OR. ( PLTRU*THRESH.LT.PL( 1 ) ) ) THEN
         RESULT( 10 ) = ULPINV
      END IF
*
      NTESTT = NTESTT + NTEST
*
*     Print out tests which fail.
*
      DO 120 J = 1, NTEST
         IF( RESULT( J ).GE.THRESH ) THEN
*
*           If this is the first test to fail,
*           print a header to the data file.
*
            IF( NERRS.EQ.0 ) THEN
               WRITE( NOUT, FMT = 9995 )'SGX'
*
*              Matrix types
*
               WRITE( NOUT, FMT = 9994 )
*
*              Tests performed
*
               WRITE( NOUT, FMT = 9992 )'orthogonal', '''', 'transpose', ( '''', I = 1, 4 )
*
            END IF
            NERRS = NERRS + 1
            IF( RESULT( J ).LT.10000.0 ) THEN
               WRITE( NOUT, FMT = 9989 )NPTKNT, MPLUSN, J, RESULT( J )
            ELSE
               WRITE( NOUT, FMT = 9988 )NPTKNT, MPLUSN, J, RESULT( J )
            END IF
         END IF
*
  120 CONTINUE
*
  130 CONTINUE
      GO TO 80
  140 CONTINUE
*
  150 CONTINUE
*
*     Summary
*
      CALL ALASVM( 'SGX', NOUT, NERRS, NTESTT, 0 )
*
      WORK( 1 ) = MAXWRK
*
      RETURN
*
 9999 FORMAT( ' SDRGSX: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ')' )
*
 9998 FORMAT( ' SDRGSX: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', Input Example #', I2, ')' )
*
 9997 FORMAT( ' SDRGSX: SGET53 returned INFO=', I1, ' for eigenvalue ',
     $      I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ')' )
*
 9996 FORMAT( ' SDRGSX: S not in Schur form at eigenvalue ', I6, '.',
     $      / 9X, 'N=', I6, ', JTYPE=', I6, ')' )
*
 9995 FORMAT( / 1X, A3, ' -- Real Expert Generalized Schur form',
     $      ' problem driver' )
*
 9994 FORMAT( 'Input Example' )
*
 9993 FORMAT( ' Matrix types: ', /
     $      '  1:  A is a block diagonal matrix of Jordan blocks ',
     $      'and B is the identity ', / '      matrix, ',
     $      / '  2:  A and B are upper triangular matrices, ',
     $      / '  3:  A and B are as type 2, but each second diagonal ',
     $      'block in A_11 and ', /
     $      '      each third diagonal block in A_22 are 2x2 blocks,',
     $      / '  4:  A and B are block diagonal matrices, ',
     $      / '  5:  (A,B) has potentially close or common ',
     $      'eigenvalues.', / )
*
 9992 FORMAT( / ' Tests performed:  (S is Schur, T is triangular, ',
     $      'Q and Z are ', A, ',', / 19X,
     $      ' a is alpha, b is beta, and ', A, ' means ', A, '.)',
     $      / '  1 = | A - Q S Z', A,
     $      ' | / ( |A| n ulp )      2 = | B - Q T Z', A,
     $      ' | / ( |B| n ulp )', / '  3 = | I - QQ', A,
     $      ' | / ( n ulp )             4 = | I - ZZ', A,
     $      ' | / ( n ulp )', / '  5 = 1/ULP  if A is not in ',
     $      'Schur form S', / '  6 = difference between (alpha,beta)',
     $      ' and diagonals of (S,T)', /
     $      '  7 = 1/ULP  if SDIM is not the correct number of ',
     $      'selected eigenvalues', /
     $      '  8 = 1/ULP  if DIFEST/DIFTRU > 10*THRESH or ',
     $      'DIFTRU/DIFEST > 10*THRESH',
     $      / '  9 = 1/ULP  if DIFEST <> 0 or DIFTRU > ULP*norm(A,B) ',
     $      'when reordering fails', /
     $      ' 10 = 1/ULP  if PLEST/PLTRU > THRESH or ',
     $      'PLTRU/PLEST > THRESH', /
     $      '    ( Test 10 is only for input examples )', / )
 9991 FORMAT( ' Matrix order=', I2, ', type=', I2, ', a=', E10.3,
     $      ', order(A_11)=', I2, ', result ', I2, ' is ', 0P, F8.2 )
 9990 FORMAT( ' Matrix order=', I2, ', type=', I2, ', a=', E10.3,
     $      ', order(A_11)=', I2, ', result ', I2, ' is ', 0P, E10.3 )
 9989 FORMAT( ' Input example #', I2, ', matrix order=', I4, ',',
     $      ' result ', I2, ' is', 0P, F8.2 )
 9988 FORMAT( ' Input example #', I2, ', matrix order=', I4, ',',
     $      ' result ', I2, ' is', 1P, E10.3 )
*
*     End of SDRGSX
*
      END
