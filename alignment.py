import numpy as np
import token
from io import BytesIO
from tokenize import tokenize

class Alignment(object):
    # Creates a simple similarity furction from a match and mismatch cost
    @staticmethod
    def simple_similarity(match=1, mismatch=-1):
        return lambda a, b: match if a == b else mismatch

    @staticmethod
    def char_aln2str(c):
        if c:
            return c
        else:
            return "-"

    def __init__(self, sequence_1, sequence_2):
        self._seq1 = sequence_1
        self._seq2 = sequence_2
        self.gap_penalty = -1
        self.method = "global"
        self._mat = None
        self._score = None
        self._seq1_aln = None
        self._seq2_aln = None
        self._similarity = Alignment.simple_similarity()
        self._aln2str = Alignment.char_aln2str
        self._alneq = lambda a, b: a == b

    def __repr__(self):
        return 'Alignment(score="{}",method="{}",gap_penalty="{}")'.format(self._score, self.method, self.gap_penalty)

    def __str__(self):
        if not self._score:
            # unaligned
            s = self.__repr__() + ", (unaligned)\n"
            s += "> {}\n> {}".format(''.join(str(c) for c in self._seq1), ''.join(str(c) for c in self._seq2))
        else:
            #TODO: adjust for different print sizes
            s = '{}, (aligned)\n> {}\n  {}\n> {}'.format(self.__repr__(),
                     ''.join(self._aln2str(c) for c in self._seq1_aln),
                     ''.join('|' if self._alneq(self._seq1_aln[i], self._seq2_aln[i]) else " " for i in range(len(self._seq1_aln))),
                     ''.join(self._aln2str(c) for c in self._seq2_aln))

        return s

    def seq1(self, seq=None):
        if seq:
            self._seq=seq
            return self
        return self._seq1

    def seq2(self, seq=None):
        if seq:
            self._seq=seq
            return self
        return self._seq1

    def alignment(self):
        return (self._seq1_aln, self._seq2_aln)

    @property
    def score(self):
        return self._score

    def similarity(self, func=None):
        if func:
            self._similarity = func
            return self
        return self._similarity


    def align(self):
        self._mat = np.zeros((len(self._seq1)+1, len(self._seq2)+1))

        _max_score = 0
        _max_i = 0
        _max_j = 0

        # initialize first row and column, if doing global
        if self.method == "global":
            for i in range(len(self._seq1)):
                self._mat[i+1][0] = self._mat[i][0] + self.gap_penalty

            for j in range(len(self._seq2)):
                self._mat[0][j+1] = self._mat[0][j] + self.gap_penalty

        # fill in the table
        for i in range(1,len(self._seq1)+1):
            for j in range(1, len(self._seq2)+1):

                left = self._mat[i][j-1] + self.gap_penalty
                up   = self._mat[i-1][j] + self.gap_penalty
                diag = self._mat[i-1][j-1] + self._similarity(self._seq1[i-1], self._seq2[j-1]) # compare current characters, remembering i and j index the table, one less than the strings

                score = max(left, up, diag)

                if self.method == "local" and score < 0:
                    score = 0

                # keep track of the maximum score position, for local alignments
                if score > _max_score:
                    _max_i = i
                    _max_j = j
                    _max_score = score

                self._mat[i][j] = score

        # Backtrack to find the 
        self._seq1_aln = []
        self._seq2_aln = []

        if self.method == "global":
            b_i = len(self._seq1)
            b_j = len(self._seq2)
        else:
            b_i = _max_i
            b_j = _max_j

        while b_i > 0 or b_j > 0:
            # check to see which box this came from
            left = self._mat[b_i][b_j-1] + self.gap_penalty
            up   = self._mat[b_i-1][b_j] + self.gap_penalty
            diag = self._mat[b_i-1][b_j-1] + self._similarity(self._seq1[b_i-1], self._seq2[b_j-1]) 

            # find out which sequnce to take elements from. None indicates a gap
            if self._mat[b_i][b_j] == diag:
                self._seq1_aln.insert(0, self._seq1[b_i-1])
                self._seq2_aln.insert(0, self._seq2[b_j-1])
                b_i -= 1
                b_j -= 1
            elif self._mat[b_i][b_j] == left:
                self._seq1_aln.insert(0, None)
                self._seq2_aln.insert(0, self._seq2[b_j-1])
                b_j -= 1
            else:
                self._seq1_aln.insert(0, self._seq1[b_i-1])
                self._seq2_aln.insert(0, None)
                b_i -= 1

            # end local alignments early
            if self.method == "local" and self._mat[b_i][b_j] == 0:
                break

        self._score = _max_score if self.method=="local" else self._mat[-1][-1]

        return self


class PySourceAlignment(Alignment):
    @staticmethod
    def _tokenize(s):
        gen = tokenize(BytesIO(s.encode('utf-8')).readline)
        # skip encoding
        next(gen)
        for tok in gen:
            yield tok

    @staticmethod
    def simple_token_similarity(match=3, mismatch=-2):
        return lambda a, b: match if a.type == b.type else mismatch

    TOK_STR_DICT = {
        token.ENDMARKER : "Z",
        token.NAME : "A",
        token.NUMBER : "1",
        token.STRING : "S",
        token.NEWLINE : "n",
        token.INDENT : "I",
        token.DEDENT : "D",
        token.LPAR : "(",
        token.RPAR : ")",
        token.LSQB : "{",
        token.RSQB : "}",
        token.COLON : ":",
        token.COMMA : ",",
        token.SEMI : ";",
        token.PLUS : "+",
        token.MINUS : "-",
        token.STAR : "*",
        token.SLASH : "/",
        token.VBAR : "|",
        token.AMPER : "&",
        token.LESS : "<",
        token.GREATER : ">",
        token.EQUAL : "=",
        token.DOT : ".",
        token.PERCENT : "%",
        token.LBRACE : "{",
        token.RBRACE : "}",
        token.EQEQUAL : "E",
        token.NOTEQUAL : "N",
        token.LESSEQUAL : "L",
        token.GREATEREQUAL : "G",
        token.TILDE : "~",
        token.CIRCUMFLEX : "^",
        token.LEFTSHIFT : "L",
        token.RIGHTSHIFT : "R",
        token.DOUBLESTAR : "s",
        token.PLUSEQUAL : "P",
        token.MINEQUAL : "M",
        token.STAREQUAL : "t",
        token.SLASHEQUAL : "u",
        token.PERCENTEQUAL : "p",
        token.AMPEREQUAL : "a",
        token.VBAREQUAL : "B",
        token.CIRCUMFLEXEQUAL : "C",
        token.LEFTSHIFTEQUAL : "h",
        token.RIGHTSHIFTEQUAL : "H",
        token.DOUBLESTAREQUAL : "T",
        token.DOUBLESLASH : "b",
        token.DOUBLESLASHEQUAL : "U",
        token.AT : "@",
        token.ATEQUAL : "j",
        token.RARROW : "W",
        token.ELLIPSIS : "e",
        token.OP : "!",
        token.AWAIT : "k",
        token.ASYNC : "K",
        token.ERRORTOKEN : "?",
        token.N_TOKENS : "x",
        token.NT_OFFSET : "X",
        # FIXME: ENCODING is not a real token
        59 : "`"

    }

    @staticmethod
    def _tok_aln2str(t):
        if t:
            return PySourceAlignment.TOK_STR_DICT[t.type]
        else:
            return " "


    def __init__(self, source1, source2):
        super().__init__(list(PySourceAlignment._tokenize(source1)), list(PySourceAlignment._tokenize(source2)))
        self._aln2str = PySourceAlignment._tok_aln2str
        self._similarity = PySourceAlignment.simple_token_similarity()
        self._alneq = lambda t1, t2 : t1 and t2 and (t1.type == t2.type)




