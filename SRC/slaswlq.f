*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim

*  Authors:
*  ========

*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.

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

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N, MB, NB, LWORK, LDT;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), WORK( * ), T( LDT, * )
      // ..

*  =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, II, KK, CTR, MINMN, LWMIN;
      // ..
      // .. EXTERNAL FUNCTIONS ..
      bool               LSAME;
      // EXTERNAL LSAME
      REAL               SROUNDUP_LWORK
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. EXTERNAL SUBROUTINES ..
      // EXTERNAL SGELQT, SGEQRT, STPLQT, STPQRT, XERBLA
      // ..
      // .. INTRINSIC FUNCTIONS ..
      // INTRINSIC MAX, MIN, MOD
      // ..
      // .. EXECUTABLE STATEMENTS ..

      // TEST THE INPUT ARGUMENTS

      INFO = 0

      LQUERY = ( LWORK == -1 )

      MINMN = MIN( M, N )
      if ( MINMN == 0 ) {
        LWMIN = 1
      } else {
        LWMIN = M*MB
      }

      if ( M.LT.0 ) {
        INFO = -1
      } else if ( N.LT.0 .OR. N.LT.M ) {
        INFO = -2
      } else if ( MB.LT.1 .OR. ( MB.GT.M .AND. M.GT.0 ) ) {
        INFO = -3
      } else if ( NB.LE.0 ) {
        INFO = -4
      } else if ( LDA.LT.MAX( 1, M ) ) {
        INFO = -6
      } else if ( LDT.LT.MB ) {
        INFO = -8
      } else if ( LWORK.LT.LWMIN .AND. (.NOT.LQUERY) ) {
        INFO = -10
      }
      if ( INFO == 0 ) {
        WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      }

      if ( INFO != 0 ) {
        xerbla('SLASWLQ', -INFO );
        RETURN
      } else if ( LQUERY ) {
        RETURN
      }

      // Quick return if possible

      if ( MINMN == 0 ) {
        RETURN
      }

      // The LQ Decomposition

      if ( (M.GE.N) .OR. (NB.LE.M) .OR. (NB.GE.N) ) {
        sgelqt(M, N, MB, A, LDA, T, LDT, WORK, INFO );
        RETURN
      }

      KK = MOD((N-M),(NB-M))
      II = N-KK+1

      // Compute the LQ factorization of the first block A(1:M,1:NB)

      sgelqt(M, NB, MB, A(1,1), LDA, T, LDT, WORK, INFO );
      CTR = 1

      DO I = NB+1, II-NB+M, (NB-M)

        // Compute the QR factorization of the current block A(1:M,I:I+NB-M)

        stplqt(M, NB-M, 0, MB, A(1,1), LDA, A( 1, I ), LDA, T(1, CTR * M + 1), LDT, WORK, INFO );
        CTR = CTR + 1
      }

      // Compute the QR factorization of the last block A(1:M,II:N)

      if ( II.LE.N ) {
        stplqt(M, KK, 0, MB, A(1,1), LDA, A( 1, II ), LDA, T(1, CTR * M + 1), LDT, WORK, INFO );
      }

      WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      RETURN

      // End of SLASWLQ

      }
