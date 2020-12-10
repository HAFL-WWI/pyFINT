######################################################################
# This script is based on the software FINT (C++ implementation v1.10 
# from July 2017; (C) ecorisQ - Luuk Dorren, Nicolas Zuanon)
#
# Author: Christoph Schaller, BFH-HAFL, December 2020
#
# Script with class for parsing DBH functions with python based syntax. 
#
# Basing on https://raw.githubusercontent.com/pyparsing/pyparsing/master/examples/fourFn.py
# and https://raw.githubusercontent.com/pyparsing/pyparsing/master/examples/SimpleCalc.py
#######################################################################
from pyparsing import Literal,Word,Group,\
    ZeroOrMore,Forward,alphas,alphanums,Regex,ParseException,\
    CaselessKeyword, Suppress
import math
import operator

class DbhParser:
    expr = None
    exprStack = []
    variables = {}
    def __init__(self):
        self.bnf = self.BNF()
    
    def clear_var(self):
        self.variables.clear()
    
    def define_var(self, var, value):
        self.variables[var]=value

    def set_expr(self, expression):
        self.expr=expression

    def eval(self):
        self.exprStack[:] = []
        
        try:
            results = self.bnf.parseString( self.expr, parseAll=True )
            val = self.evaluateStack( self.exprStack[:] )
        except ParseException as pe:
            print("failed parse:", str(pe))
        except Exception as e:
            print("failed eval:", str(e))
        else:
            return val

    def evaluateStack(self, s ):
        op = s.pop()
        if op == 'unary -':
            return -self.evaluateStack( s )
        if op in "+-*/^":
            op2 = self.evaluateStack( s )
            op1 = self.evaluateStack( s )
            return self.opn[op]( op1, op2 )
        elif op == "PI":
            return math.pi # 3.1415926535
        elif op == "E":
            return math.e  # 2.718281828
        elif op in self.fn:
            return self.fn[op]( self.evaluateStack( s ) )
        elif op[0].isalpha():
            if op in self.variables:
                return self.variables[op]
            raise Exception("invalid identifier '%s'" % op)
        else:
            return float( op )

    def pushFirst(self, strg, loc, toks ):
        self.exprStack.append( toks[0] )
    def pushUMinus(self, strg, loc, toks ):
        for t in toks:
            if t == '-':
                self.exprStack.append( 'unary -' )
                #~ exprStack.append( '-1' )
                #~ exprStack.append( '*' )
            else:
                break

    bnf = None
    def BNF(self):
        """
        expop   :: '^'
        multop  :: '*' | '/'
        addop   :: '+' | '-'
        integer :: ['+' | '-'] '0'..'9'+
        atom    :: PI | E | real | fn '(' expr ')' | '(' expr ')'
        factor  :: atom [ expop factor ]*
        term    :: factor [ multop factor ]*
        expr    :: term [ addop term ]*
        """
        if not self.bnf:
            point = Literal( "." )
            # use CaselessKeyword for e and pi, to avoid accidentally matching
            # functions that start with 'e' or 'pi' (such as 'exp'); Keyword
            # and CaselessKeyword only match whole words
            e     = CaselessKeyword( "E" )
            pi    = CaselessKeyword( "PI" )
            #~ fnumber = Combine( Word( "+-"+nums, nums ) +
                            #~ Optional( point + Optional( Word( nums ) ) ) +
                            #~ Optional( e + Word( "+-"+nums, nums ) ) )
            fnumber = Regex(r"[+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?")
            ident = Word(alphas, alphanums+"_$")

            plus, minus, mult, div = map(Literal, "+-*/")
            lpar, rpar = map(Suppress, "()")
            addop  = plus | minus
            multop = mult | div
            expop = Literal( "^" )

            expr = Forward()
            atom = ((0,None)*minus + ( pi | e | fnumber | ident + lpar + expr + rpar | ident ).setParseAction( self.pushFirst ) |
                    Group( lpar + expr + rpar )).setParseAction(self.pushUMinus)

            # by defining exponentiation as "atom [ ^ factor ]..." instead of "atom [ ^ atom ]...", we get right-to-left exponents, instead of left-to-righ
            # that is, 2^3^2 = 2^(3^2), not (2^3)^2.
            factor = Forward()
            factor << atom + ZeroOrMore( ( expop + factor ).setParseAction( self.pushFirst ) )

            term = factor + ZeroOrMore( ( multop + factor ).setParseAction( self.pushFirst ) )
            expr << term + ZeroOrMore( ( addop + term ).setParseAction( self.pushFirst ) )
            self.bnf = expr
        return self.bnf

    # map operator symbols to corresponding arithmetic operations
    epsilon = 1e-12
    opn = { "+" : operator.add,
            "-" : operator.sub,
            "*" : operator.mul,
            "/" : operator.truediv,
            "^" : operator.pow }
    fn  = { "sin" : math.sin,
            "cos" : math.cos,
            "tan" : math.tan,
            "exp" : math.exp,
            "ceil": math.ceil,
            "floor": math.floor,
            "ln": math.log,
            "abs" : abs,
            "trunc" : lambda a: int(a),
            "round" : round,
            "sgn" : lambda a: (a > epsilon) - (a < -epsilon) }
