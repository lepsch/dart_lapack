*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*> Short-Wide LQ (SWLQ) performs LQ by a sequence of orthogonal transformations,
*> representing Q as a product of other orthogonal matrices
*>   Q = Q(1) * Q(2) * . . . * Q(k)
*> where each Q(i) zeros out upper diagonal entries of a block of NB rows of A:
*>   Q(1) zeros out the upper diagonal entries of rows 1:NB of A
*>   Q(2) zeros out the bottom MB-N rows of rows [1:M,NB+1:2*NB-M] of A
*>   Q(3) zeros out the bottom MB-N rows of rows [1:M,2*NB-M+1:3*NB-2*M] of A
*>   . . .
*>
*> Q(1) is computed by GELQT, which represents Q(1) by Householder vectors
*> stored under the diagonal of rows 1:MB of A, and by upper triangular
*> block reflectors, stored in array T(1:LDT,1:N).
*> For more information see Further Details in GELQT.
*>
*> Q(i) for i>1 is computed by TPLQT, which represents Q(i) by Householder vectors
*> stored in columns [(i-1)*(NB-M)+M+1:i*(NB-M)+M] of A, and by upper triangular
*> block reflectors, stored in array T(1:LDT,(i-1)*M+1:i*M).
*> The last Q(k) may use fewer rows.
*> For more information see Further Details in TPQRT.
*>
*> For more details of the overall algorithm, see the description of
*> Sequential TSQR in Section 2.2 of [1].
*>
*> [1] “Communication-Optimal Parallel and Sequential QR and LU Factorizations,”
*>     J. Demmel, L. Grigori, M. Hoemmen, J. Langou,
*>     SIAM J. Sci. Comput, vol. 34, no. 1, 2012
*> \endverbatim
*>
*> \ingroup laswlq
*>
*  =====================================================================
      SUBROUTINE SLASWLQ( M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
*
*     .. Scalar Arguments ..
      int                INFO, LDA, M, N, MB, NB, LWORK, LDT;
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), WORK( * ), T( LDT, * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      bool               LQUERY;
      int                I, II, KK, CTR, MINMN, LWMIN;
*     ..
*     .. EXTERNAL FUNCTIONS ..
      bool               LSAME;
      // EXTERNAL LSAME
      REAL               SROUNDUP_LWORK
      // EXTERNAL SROUNDUP_LWORK
*     ..
*     .. EXTERNAL SUBROUTINES ..
      // EXTERNAL SGELQT, SGEQRT, STPLQT, STPQRT, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      // INTRINSIC MAX, MIN, MOD
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT ARGUMENTS
*
      INFO = 0
*
      LQUERY = ( LWORK.EQ.-1 )
*
      MINMN = MIN( M, N )
      IF( MINMN.EQ.0 ) THEN
        LWMIN = 1
      ELSE
        LWMIN = M*MB
      END IF
*
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 .OR. N.LT.M ) THEN
        INFO = -2
      ELSE IF( MB.LT.1 .OR. ( MB.GT.M .AND. M.GT.0 ) ) THEN
        INFO = -3
      ELSE IF( NB.LE.0 ) THEN
        INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -6
      ELSE IF( LDT.LT.MB ) THEN
        INFO = -8
      ELSE IF( LWORK.LT.LWMIN .AND. (.NOT.LQUERY) ) THEN
        INFO = -10
      END IF
      IF( INFO.EQ.0 ) THEN
        WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      END IF
*
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'SLASWLQ', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
*
*     Quick return if possible
*
      IF( MINMN.EQ.0 ) THEN
        RETURN
      END IF
*
*     The LQ Decomposition
*
      IF( (M.GE.N) .OR. (NB.LE.M) .OR. (NB.GE.N) ) THEN
        CALL SGELQT( M, N, MB, A, LDA, T, LDT, WORK, INFO )
        RETURN
      END IF
*
      KK = MOD((N-M),(NB-M))
      II = N-KK+1
*
*     Compute the LQ factorization of the first block A(1:M,1:NB)
*
      CALL SGELQT( M, NB, MB, A(1,1), LDA, T, LDT, WORK, INFO )
      CTR = 1
*
      DO I = NB+1, II-NB+M, (NB-M)
*
*       Compute the QR factorization of the current block A(1:M,I:I+NB-M)
*
        CALL STPLQT( M, NB-M, 0, MB, A(1,1), LDA, A( 1, I ), LDA, T(1, CTR * M + 1), LDT, WORK, INFO )
        CTR = CTR + 1
      END DO
*
*     Compute the QR factorization of the last block A(1:M,II:N)
*
      IF( II.LE.N ) THEN
        CALL STPLQT( M, KK, 0, MB, A(1,1), LDA, A( 1, II ), LDA, T(1, CTR * M + 1), LDT, WORK, INFO )
      END IF
*
      WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      RETURN
*
*     End of SLASWLQ
*
      END
