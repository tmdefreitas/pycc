import numpy as np

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
                     ''.join('|' if (self._seq1_aln[i] == self._seq2_aln[i]) else " " for i in range(len(self._seq1_aln))),
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


