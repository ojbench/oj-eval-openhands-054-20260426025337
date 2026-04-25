
#include <bits/stdc++.h>
using namespace std;

class Term {
public:
    int a, b, c, d;  // coefficient, x^b, sin^c(x), cos^d(x)
    
    Term(int a = 0, int b = 0, int c = 0, int d = 0) : a(a), b(b), c(c), d(d) {}
    
    bool operator==(const Term& other) const {
        return b == other.b && c == other.c && d == other.d;
    }
    
    bool operator!=(const Term& other) const {
        return !(*this == other);
    }
    
    bool operator<(const Term& other) const {
        if (b != other.b) return b > other.b;
        if (c != other.c) return c > other.c;
        return d > other.d;
    }
};

class Poly {
public:
    vector<Term> terms;
    
    Poly() {}
    Poly(const vector<Term>& t) : terms(t) {}
    
    void simplify() {
        if (terms.empty()) return;
        
        sort(terms.begin(), terms.end());
        vector<Term> simplified;
        
        for (const Term& t : terms) {
            if (t.a == 0) continue;
            
            if (!simplified.empty() && simplified.back() == t) {
                simplified.back().a += t.a;
            } else {
                simplified.push_back(t);
            }
        }
        
        terms.clear();
        for (const Term& t : simplified) {
            if (t.a != 0) {
                terms.push_back(t);
            }
        }
    }
    
    Poly operator+(const Poly& other) const {
        Poly result;
        result.terms = terms;
        result.terms.insert(result.terms.end(), other.terms.begin(), other.terms.end());
        result.simplify();
        return result;
    }
    
    Poly operator-(const Poly& other) const {
        Poly result;
        result.terms = terms;
        for (Term t : other.terms) {
            t.a = -t.a;
            result.terms.push_back(t);
        }
        result.simplify();
        return result;
    }
    
    Poly operator*(const Poly& other) const {
        Poly result;
        for (const Term& t1 : terms) {
            for (const Term& t2 : other.terms) {
                Term new_term(t1.a * t2.a, t1.b + t2.b, t1.c + t2.c, t1.d + t2.d);
                result.terms.push_back(new_term);
            }
        }
        result.simplify();
        return result;
    }
    
    Poly copy() const {
        Poly result;
        result.terms = terms;
        return result;
    }
    
    Poly derivate() const {
        Poly result;
        for (const Term& t : terms) {
            if (t.a == 0) continue;
            
            // Derivative of a * x^b * sin^c(x) * cos^d(x)
            // Using product rule: (uvw)' = u'vw + uv'w + uvw'
            
            // u = a * x^b, v = sin^c(x), w = cos^d(x)
            
            // u' * v * w
            if (t.b > 0) {
                Term term1(t.a * t.b, t.b - 1, t.c, t.d);
                result.terms.push_back(term1);
            }
            
            // u * v' * w
            if (t.c > 0) {
                // (sin^c(x))' = c * sin^(c-1)(x) * cos(x)
                Term term2(t.a * t.c, t.b, t.c - 1, t.d + 1);
                result.terms.push_back(term2);
            }
            
            // u * v * w'
            if (t.d > 0) {
                // (cos^d(x))' = -d * sin(x) * cos^(d-1)(x)
                Term term3(-t.a * t.d, t.b, t.c + 1, t.d - 1);
                result.terms.push_back(term3);
            }
        }
        result.simplify();
        return result;
    }
};

class Frac {
public:
    Poly p, q;
    
    Frac() {}
    Frac(int x) {
        p.terms.push_back(Term(x, 0, 0, 0));
        q.terms.push_back(Term(1, 0, 0, 0));
    }
    Frac(const Term& t) {
        p.terms.push_back(t);
        q.terms.push_back(Term(1, 0, 0, 0));
    }
    Frac(const Poly& _p, const Poly& _q) : p(_p), q(_q) {}
    
    Frac operator+(const Frac& other) const {
        Frac result;
        result.p = p * other.q + other.p * q;
        result.q = q * other.q;
        return result;
    }
    
    Frac operator-(const Frac& other) const {
        Frac result;
        result.p = p * other.q - other.p * q;
        result.q = q * other.q;
        return result;
    }
    
    Frac operator*(const Frac& other) const {
        Frac result;
        result.p = p * other.p;
        result.q = q * other.q;
        return result;
    }
    
    Frac operator/(const Frac& other) const {
        Frac result;
        result.p = p * other.q;
        result.q = q * other.p;
        return result;
    }
    
    Frac derivate() const {
        Frac result;
        result.p = p.derivate() * q - q.derivate() * p;
        result.q = q * q;
        return result;
    }
    
    string toString() const {
        string p_str = polyToString(p);
        string q_str = polyToString(q);
        
        if (p_str == "0") return "0";
        if (q_str == "1") return p_str;
        
        bool p_need_paren = p.terms.size() > 1;
        bool q_need_paren = q.terms.size() > 1;
        
        string result;
        if (p_need_paren) result += "(";
        result += p_str;
        if (p_need_paren) result += ")";
        result += "/";
        if (q_need_paren) result += "(";
        result += q_str;
        if (q_need_paren) result += ")";
        
        return result;
    }
    
private:
    string polyToString(const Poly& poly) const {
        if (poly.terms.empty()) return "0";
        
        string result;
        for (size_t i = 0; i < poly.terms.size(); i++) {
            const Term& t = poly.terms[i];
            
            if (i == 0) {
                if (t.a < 0) result += "-";
            } else {
                if (t.a < 0) result += "-";
                else result += "+";
            }
            
            int coeff = abs(t.a);
            if (coeff != 1 || (t.b == 0 && t.c == 0 && t.d == 0)) {
                result += to_string(coeff);
            }
            
            if (t.b > 0) {
                result += "x";
                if (t.b > 1) result += "^" + to_string(t.b);
            }
            
            if (t.c > 0) {
                result += "sin";
                if (t.c > 1) result += "^" + to_string(t.c);
                result += "x";
            }
            
            if (t.d > 0) {
                result += "cos";
                if (t.d > 1) result += "^" + to_string(t.d);
                result += "x";
            }
        }
        
        return result;
    }
};

int get_number(const string& s, int l, int r) {
    bool negative = false;
    int num = 0;
    bool has_digit = false;
    
    for (int i = l; i < r; i++) {
        if (s[i] == '-') {
            negative = true;
        } else if (isdigit(s[i])) {
            num = num * 10 + (s[i] - '0');
            has_digit = true;
        }
    }
    
    if (!has_digit) {
        return negative ? -1 : 1;
    }
    
    return negative ? -num : num;
}

Term get_term(const string& s, int l, int r) {
    int coeff = 1;
    int b = 0, c = 0, d = 0;
    
    // Extract coefficient
    int i = l;
    bool negative = false;
    if (i < r && s[i] == '-') {
        negative = true;
        i++;
    }
    
    // Parse number part
    int num = 0;
    bool has_digit = false;
    while (i < r && isdigit(s[i])) {
        num = num * 10 + (s[i] - '0');
        i++;
        has_digit = true;
    }
    
    if (has_digit) {
        coeff = num;
    } else {
        coeff = 1;
    }
    
    if (negative) coeff = -coeff;
    
    // Parse variables and functions
    while (i < r) {
        if (s[i] == 'x' && (i == l || (s[i-1] != 'n' && s[i-1] != 's'))) {
            // x power
            if (i + 1 < r && s[i+1] == '^') {
                int j = i + 2;
                int power = 0;
                while (j < r && isdigit(s[j])) {
                    power = power * 10 + (s[j] - '0');
                    j++;
                }
                b = power;
                i = j;
            } else {
                b = 1;
                i++;
            }
        } else if (i + 2 < r && s[i] == 's' && s[i+1] == 'i' && s[i+2] == 'n') {
            // sin power - format is sin^nx
            i += 3;
            if (i < r && s[i] == '^') {
                int j = i + 1;
                int power = 0;
                while (j < r && isdigit(s[j])) {
                    power = power * 10 + (s[j] - '0');
                    j++;
                }
                c = power;
                i = j;  // Skip the 'x' at the end
                if (i < r && s[i] == 'x') i++;
            } else {
                c = 1;
                if (i < r && s[i] == 'x') i++;
            }
        } else if (i + 2 < r && s[i] == 'c' && s[i+1] == 'o' && s[i+2] == 's') {
            // cos power - format is cos^nx
            i += 3;
            if (i < r && s[i] == '^') {
                int j = i + 1;
                int power = 0;
                while (j < r && isdigit(s[j])) {
                    power = power * 10 + (s[j] - '0');
                    j++;
                }
                d = power;
                i = j;  // Skip the 'x' at the end
                if (i < r && s[i] == 'x') i++;
            } else {
                d = 1;
                if (i < r && s[i] == 'x') i++;
            }
        } else {
            i++;
        }
    }
    
    return Term(coeff, b, c, d);
}

Frac dfs(const string& s, int l, int r) {
    // Handle parentheses - find rightmost + or - for left-to-right evaluation
    int balance = 0;
    int last_op = -1;
    for (int i = l; i < r; i++) {
        if (s[i] == '(') balance++;
        else if (s[i] == ')') balance--;
        else if (balance == 0 && (s[i] == '+' || s[i] == '-')) {
            if (i == l) continue; // Skip leading negative
            last_op = i;
        }
    }
    
    if (last_op != -1) {
        Frac left = dfs(s, l, last_op);
        Frac right = dfs(s, last_op + 1, r);
        if (s[last_op] == '+') return left + right;
        else return left - right;
    }
    
    // Handle multiplication and division
    balance = 0;
    for (int i = l; i < r; i++) {
        if (s[i] == '(') balance++;
        else if (s[i] == ')') balance--;
        else if (balance == 0 && (s[i] == '*' || s[i] == '/')) {
            Frac left = dfs(s, l, i);
            Frac right = dfs(s, i + 1, r);
            if (s[i] == '*') return left * right;
            else return left / right;
        }
    }
    
    // Handle parentheses
    if (s[l] == '(' && s[r-1] == ')') {
        return dfs(s, l + 1, r - 1);
    }
    
    // Single term
    Term t = get_term(s, l, r);
    return Frac(t);
}

void solve(const string& s, int n) {
    Frac f = dfs(s, 0, n);
    Frac f_prime = f.derivate();
    
    cout << f.toString() << endl;
    cout << f_prime.toString() << endl;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    string expr;
    getline(cin, expr);
    
    solve(expr, expr.length());
    
    return 0;
}
