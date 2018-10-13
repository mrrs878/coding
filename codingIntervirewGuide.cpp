#include <iostream>
#include <string>
#include <stack>
#include <queue>
#include <set>
#include <list>
#include <map>
#include <unordered_map>
#include <algorithm>

using namespace std;

#define MIN(x, y)      ((x) < (y) ? (x) : (y))
#define MAX(x, y)      ((x) > (y) ? (x) : (y))
#define ALIGN(x, c);      {for(int i = 0; i<x; ++i) cout << c;}
int intlen(int _x) {    
	int len = 0;        
	if (_x == 0)        
		return 1;       
	while (_x != 0) {   
		_x /= 10;       
		++len;          
	}                   
	return len;         
}                       
string int2str(int _x) {
	string tmp, res;    
	do {                
		tmp.push_back((_x % 10) + '0');
		_x /= 10;       
	} while (_x != 0);  
	for (size_t index = tmp.size() - 1; index != -1; --index)
		res.push_back(tmp[index]);
	return res;         
}

void stackPrint(stack<int> &_tmp) {
	while (!_tmp.empty()) {
		cout << _tmp.top() << ' ';
		_tmp.pop();
	}
}

template<class type>
void vectorPrint(vector<type> &_tmp) {
	for (vector<type>::iterator iter = _tmp.begin(); iter != _tmp.end(); ++iter)
		cout << *iter << ' ';
	cout << endl;
}

void listPrint(list<int> &_tmp) {
	for (list<int>::iterator iter = _tmp.begin(); iter != _tmp.end(); ++iter)
		cout << (*iter) << ' ';
	cout << endl;
}

vector<int> generateRandomArray(int _maxsize, int _maxvalue) {
	vector<int> tmp(rand() % _maxsize);
	for (int index = 0; index < tmp.size(); ++index)
		tmp[index] = (rand() % _maxvalue);

	return tmp;
}

struct nod {
	int value;
	nod *lchild;
	nod *rchild;
	int pos;
	nod(int _value = 0, int _pos = 0, nod *_lchild = nullptr, nod *_rchild = nullptr) :
		value(_value), lchild(_lchild), rchild(_rchild), pos(_pos) {}
};
//#0  template
class type {
/****************************test code****************************/
/*

****debug:
*/
/*****************************end test****************************/
};
//#1
class stackWithGetMin {
private:
	stack<int> src;
	stack<int> help;
public:
	void POP(void) {
		if (src.top() == help.top())
			help.pop();
		src.pop();
	}
	void PUSH(int _value) {
		src.push(_value);
		if (help.empty() == true)
			help.push(_value);
		else if(_value <= help.top())
			help.push(_value);
	}
	int getMin(void) {
		return help.top();
	}
	int getTop(void) {
		return src.top();
	}
	bool isEmpty(void) {
		return src.empty();
	}
/****************************test code****************************/
/*
    stackWithGetMin s;

	s.PUSH(3); s.PUSH(4); s.PUSH(5);
	cout << "Min: " << s.getMin() << endl;
	s.PUSH(1); s.PUSH(2);
	cout << "Min: " << s.getMin() << endl;
	s.POP(); s.POP();
	cout << "Min: " << s.getMin() << endl;
*/
/****************************test code****************************/
};
//#2
class queueMadeBy2Stack {
private:
	stack<int> push;
	stack<int> pop;
public:
	void PUSH(int _value) {
		push.push(_value);
	}
	void POP(void) {
		if (pop.empty() == true) {
			while (!push.empty()) {
				pop.push(push.top());
				push.pop();
			}
		}
		pop.pop();
	}
	int getTop(void) {
		return pop.top();
	}
/****************************test code****************************/
/*	queueMadeBy2Stack s;

	s.PUSH(1); s.PUSH(2); s.PUSH(3); s.PUSH(4);
	s.POP();
	cout << s.getTop();
	s.POP();
	cout << s.getTop();
	s.PUSH(5);
	s.PUSH(6);
	s.POP();
	cout << s.getTop();
	s.POP();
	s.POP();
	cout << s.getTop();
*/
//debug: 2346
/****************************test code****************************/
};
//#3
class reverseStackByRecurvise {
private:
	stack<int> tmp;
	int getAndRmvBtm(stack<int> &_tmp) {
		int result = _tmp.top();
		_tmp.pop();
		if (_tmp.empty())
			return result;
		else {
			int btm = getAndRmvBtm(_tmp);
			tmp.push(result);
			return btm;
		}
	}
	void reverseStack(stack<int> &_tmp) {
		if (_tmp.empty())
			return;
		int btm = getAndRmvBtm(_tmp);
		reverseStack(_tmp);
		_tmp.push(btm);
	}
public:
	int getTop(void) {
		return tmp.top();
	}
	void PUSH(int _value) {
		tmp.push(_value);
	}
	void POP(void) {
		tmp.pop();
	}
	void reverseStack(void) {
		reverseStack(tmp);
	}
/****************************test code****************************/
/*	
    reverseStackByRecurvise t;

	t.PUSH(1); t.PUSH(2); t.PUSH(3); t.PUSH(4);
	t.reverseStack();
*/
/****************************test code****************************/
};
//#4  continue...
class catAndDogQueue {
	
};
//#5
class orderStackByStack {
private:
	stack<int> src;
	stack<int> tmp;
public:
	stack<int>& PRINT(void) {
		return src;
	}
	void PUSH(int _value) {
		src.push(_value);
	}
	void orderStack(void) {
		int value = src.top();

		tmp.push(value);
		src.pop();
		while (!src.empty()) {
			value = src.top();
			src.pop();
			while (!tmp.empty() && value > tmp.top()){
				src.push(tmp.top());
				tmp.pop();
			}
			tmp.push(value);
		}
		while (!tmp.empty()) {
			src.push(tmp.top());
			tmp.pop();
		}
	}
/****************************test code****************************/
/*	orderStackByStack s;

	s.PUSH(1); s.PUSH(3); s.PUSH(4); s.PUSH(5); s.PUSH(2);
	s.orderStack();
	stackPrint(s.PRINT());
*/
//debug: 54321
/****************************test code****************************/
};
//#6
class islandProblem {
private:
	int (*src)[4];
	int row, col;
	int sumNum, areaXNum, areaYNum;
	int col_start, col_end;
	char flag = 'A';
/*****union/find set*****/
	unordered_map<char, char> fatherMap;  //child -- father
	unordered_map<char, int> sizeMap;
	void makeSet(list<char> &_nodes) {
		fatherMap.clear();
		sizeMap.clear();
		for (list<char>::iterator iter = _nodes.begin(); iter != _nodes.end(); ++iter) {
			fatherMap.insert(pair<char, char>(*iter, *iter));
			sizeMap.insert(pair<char, int>(*iter, 1));
		}
	}
	char findHead(const char &_node) {
		char father = fatherMap.find(_node)->second;
		if (father != _node)
			father = findHead(father);
		fatherMap.insert(pair<char, char>(_node, father));
		return father;
	}
	bool isSameSet(const char &_x, const char &_y) {
		return (findHead(_x) == findHead(_y));
	}
	void unionSet(const char &_x, const char &_y) {
		char head_x = findHead(_x);
		char head_y = findHead(_y);
		if (head_x == head_y)
			return;
		int size_x = sizeMap.find(head_x)->second;
		int size_y = sizeMap.find(head_y)->second;
		if (size_x < size_y) {
			fatherMap.insert(pair<char, char>(head_x, head_y));
			sizeMap.insert(pair<char, int>(head_y, size_x));
		}
		else {
			fatherMap.insert(pair<char, char>(head_y, head_x));
			sizeMap.insert(pair<char, int>(head_x, size_y));
		}
	}
/**end of union/find set**/
private:
	void inFection(int _row, int _col, int(*_x)[4]) {
		if ((_row < 0) || (_row > row - 1) || (_col < col_start) || (_col > col_end - 1) 
			|| (_x[_row][_col] != 1))
			return;
		src[_row][_col] = flag;
		inFection(_row - 1, _col, _x);
		inFection( _row + 1, _col, _x);
		inFection(_row, _col - 1, _x);
		inFection(_row, _col + 1, _x);
	}
	void unionAreas(void) {
		sumNum = areaXNum + areaYNum;
		int col_lchild = (col >> 1) - 1;
		int col_rchild = (col >> 1);
		list<char> tmp;
		
		/*make union/find set*/
		for (int i = 0; i < row; ++i) {
			if (src[i][col_lchild] != 0)
				tmp.push_back(src[i][col_lchild]);
			if (src[i][col_rchild] != 0)
				tmp.push_back(src[i][col_rchild]);
		}
		makeSet(tmp);
		/*union same node*/
		for (int i = 0; i < row; ++i) {
			if ((src[i][col_lchild] != 0) && (src[i][col_rchild] != 0)) {
				if (!isSameSet(src[i][col_lchild], src[i][col_rchild]))
					unionSet(src[i][col_lchild], src[i][col_rchild]);
				--sumNum;
			}	
		}
	}
public:
	islandProblem(int (*_src)[4] = nullptr, int _row = 0, int _col = 0)
		: src(_src), row(_row), col(_col), areaXNum(0), areaYNum(0) {}
	void getNums(int row_start, int row_end, int col_start, int col_end) {
		for (int i = row_start; i < row_end; ++i) {
			for (int j = col_start; j < col_end; ++j) {
				if (src[i][j] == 1) {
					col_start == 0 ? ++areaXNum : ++areaYNum;
					inFection(i, j, src);
					flag++;
				}
			}
		}
	}
	void normalSolution(void){
		col_start = 0; col_end = col;
		getNums(0, row, 0, col);
		cout << "this array have " << areaXNum << " islands!\n";
	}
	void multiProcossorSolution(void) {
		col_start = 0; col_end = col >> 1;
		getNums(0, row, 0, col >> 1);
		col_start = col >> 1; col_end = col;
		getNums(0, row, col - (col >> 1), col);
		unionAreas();
		cout << "this array have " << sumNum << " islands!\n";
	}
/****************************test code****************************/
/*	
int x[][4] = {  {0, 1, 1, 0},
				{1, 0, 1, 0},
				{0, 1, 1, 1},
				{0 ,1, 0, 0},
				{0 ,0, 1, 1},
				{0 ,1, 0, 0},};
	islandProblem a(x, 6, 4);


	//cout << "normalSoulution: ";
	//a.normalSolution();

	cout << "multiProcossorSolution: ";
	a.multiProcossorSolution();
*/
//debug: multiProcossorSolution: this array have 4 islands!
/****************************test code****************************/
};
//#7
class trieTree {
private:
	struct trieNode {
		int path;
		int end;
		trieNode *next[26];
		trieNode(int _path = 0, int _end = 0) : path(_path), end(_end) { memset(next, 0, sizeof(next)); }
	};
	trieNode src_root;
private:
	void clearNode(trieNode *_tmp) {
		for (int i = 0; i < 26; ++i) {
			if (_tmp == nullptr)
				return;
			if (_tmp->next[i] == nullptr)
				continue;
			clearNode(_tmp->next[i]);
			delete _tmp->next[i];
		}
	}
public:
	trieTree() {}
	~trieTree() { clearNode(&src_root); }
	void insertStr(const string &_str) {
		size_t i = 0;
		int index = 0;
		trieNode *tmp = &src_root;

		for (i = 0; i < _str.length(); ++i) {
			index = _str[i] - 'a';
			if (tmp->next[index] == nullptr)
				tmp->next[index] = new trieNode();
			tmp = tmp->next[index];
			tmp->path++;
		}
		tmp->end++;
	}
	bool searchStr(const string &_str) {
		size_t i = 0;
		int index = 0;
		trieNode *tmp = &src_root;
	    
		for (i = 0; i < _str.length(); ++i) {
			index = _str[i] - 'a';
			if (tmp->next[index]) {
				tmp = tmp->next[index];
				continue;
			}
			return false;
		}
		return (tmp->end ? true : false);
	}
	void deleteStr(const string &_str) {
		size_t i = 0;
		int index = 0;
		trieNode *tmp = &src_root;

		if (searchStr(_str) == false)
			return;
		for (i = 0; i < _str.length(); ++i) {
			index = _str[i] - 'a';
			if (tmp->next[index]->path-- == 0) {
				delete tmp->next[index];
				break;
			}
			tmp = tmp->next[index];
		}
		tmp->end--;
	}
	int prefFixNumber(const string &_str) {
		size_t i = 0;
		int index = 0;
		trieNode *tmp = &src_root;

		for (i = 0; i < _str.length(); ++i) {
			index = _str[i] - 'a';
			if (tmp->next[index])
				tmp = tmp->next[index];
			else
				return 0;
		}
		return tmp->path;
	}
/****************************test code****************************/
/*
	trieTree t;

	t.insertStr("abc");
	t.insertStr("abcd");
	t.deleteStr("abc");
	cout << t.searchStr("abc") << endl;
	cout << t.prefFixNumber("ab") << endl;
	debug: 0
	       1
*/
/*****************************end test****************************/
};
//#8
class minPath {
private:
	int (*src)[4];
	int row, col;
	int path;
private:
	int getMinPathRec(int(*_src)[4], int _row, int _col) {
		int tmp = _src[_row][_col];

		if(_row == row - 1 && _col == col - 1)
			return tmp;
		if (_row == row - 1)
			return tmp + getMinPathRec(_src, _row, _col + 1);
		if (_col == col - 1)
			return tmp + getMinPathRec(_src, _row + 1, _col);
		int rchild = getMinPathRec(_src, _row, _col + 1);
		int down  = getMinPathRec(_src, _row + 1, _col);
		return ((rchild < down ? rchild : down) + tmp);
	}
	int getMinPathDP(int(*_src)[4]) {
		if (_src == nullptr)
			return 0;
		int tmp = 0;
		int (*dp)[4] = new int[row][4];
		memset(dp, 0, sizeof(dp) * row * col);
		dp[0][0] = _src[0][0];
		for (int i = 1; i < col; ++i)  //init 1st row
			dp[0][i] = dp[0][i - 1] + _src[0][i];
		for (int i = 1; i < row; ++i)  //init 1st column
			dp[i][0] = dp[i - 1][0] + _src[i][0];
		for (int i = 1; i < row; ++i) {
			for (int j = 1; j < col; ++j)
				dp[i][j] = _src[i][j] + MIN(dp[i][j - 1], dp[i - 1][j]);
		}
		tmp = dp[row - 1][col - 1];
		delete[] dp;
		return tmp;
	}
	int getMinPathDPComp(int(*_src)[4]) {
		if (_src == nullptr)
			return 0;
		int *dp = new int[col];
		memset(dp, 0, sizeof(int) * col);
		dp[0] = _src[0][0];
		for (int i = 1; i < col; ++i)
			dp[i] = _src[0][i] + dp[i - 1];
		for (int i = 1; i < row; ++i) {
			dp[0] = _src[i][0] + dp[0];
			for (int j = 1; j < col; ++j) 
				dp[j] = _src[i][j] + MIN(dp[j], dp[j - 1]);
		}
		int tmp = dp[col - 1];
		delete dp;
		return tmp;
	}
public:
	minPath(int(*_src)[4] = nullptr, int _row = 0, int _col = 0) 
		: src(_src), row(_row), col(_col), path(0) {}
	int getMinPathRec(void) {
		return getMinPathRec(src, 0, 0);
	}
	int getMinPathDP(void) {
		return getMinPathDP(src);
	}
	int getMinPathDPComp(void) {
		return getMinPathDPComp(src);
	}
/****************************test code****************************/
/*
    int src[][4] = { {1, 3, 5, 9},
	                 {8, 1, 3, 4},
					 {5, 0, 6, 1},
					 {8, 8, 4, 0}, };
	minPath m(src, 4, 4);

	cout << "getMinPathRec: " << m.getMinPathRec() << endl;
	cout << "getMinPathDP: " << m.getMinPathDP() << endl;
	cout << "getMinPathDPComp: " << m.getMinPathDPComp() << endl;
****debug: getMinPathRec: 12
           getMinPathDP: 12
		   getMinPathDPComp: 12
*/
/*****************************end test****************************/
};
//#9
class isSum {
private:
	int *arr;
	int size;
	int aim;
private:
	bool judgeIsSumRec(int *_arr, int _index, int _sum, int _aim) {
		if (_index == size)
			return _sum == _aim;
		return judgeIsSumRec(_arr, _index + 1, _sum + _arr[_index], _aim) || 
			   judgeIsSumRec(_arr, _index + 1, _sum               , _aim);
	}
	bool judgeIsSumDP(int *_arr) {
		int sum = 0;  //aim = [0, sum]
		for (int i = 0; i < size; ++i)
			sum += _arr[i];
		if (sum < aim)
			return false;

		bool tmp = false;
		bool *dp = new bool[sum + 1];
		memset(dp, false, sizeof(bool) * (sum + 1));
		dp[aim] = true;
		for (int i = size - 1; i >= 0; --i) {
			for (int j = 0; j < sum; ++j) {
				if (j + _arr[i] < sum)
					dp[j] = dp[j + 1] || dp[j + _arr[i]];
				else 
					dp[j] = (j == aim ? true : false);
				tmp |= dp[j];
			}
		}
		delete[] dp;
		return tmp;
	}
public:
	isSum(int *_arr = nullptr, int _size = 0, int _aim = 0) : arr(_arr), aim(_aim), size(_size) {}
	bool judgeIsSumRec(void) {
		return judgeIsSumRec(arr, 0, 0, aim);
	}
	bool judgeIsSumDP(void) {
		return judgeIsSumDP(arr);
	}
/****************************test code****************************/
/*
	int src[4] = { 1, 3, 5, 7 };
	int aim = 15;
	isSum a(src, 4, aim);

	cout << "judgeIsSum " << aim << (a.judgeIsSum() ? ": true" : ": false") << endl;
	aim = 18;
	isSum b(src, 4, aim);
	cout << "judgeIsSum " << aim << (b.judgeIsSum() ? ": true" : ": false") << endl;
	cout << "judgeIsSumDP " << aim << (b.judgeIsSumDP()  ? ": true" : ": false") << endl
****debug: judgeIsSum 18 : false
		   judgeIsSum 7 : true
		   judgeIsSumDP: 7: true
*/
/*****************************end test****************************/
};
//#10
class fibonacciList {
private:
	int n;
private:
	int getFibonacciNN(int _n) {
		if (_n == 1 || _n == 2)
			return 1;
		return getFibonacciNN(_n - 1) + getFibonacciNN(_n - 2);
	}
public:
	fibonacciList(int _n = 0) : n(_n) {}
	int getFibonacciNN(void) {
		return getFibonacciNN(n);
	}
	int getFibonacciN(void) {
		int res = 1, tmp = 0, pre = 1;

		if (n == 1 || n == 2)
			return 1;
		for (int i = 3; i <= n; ++i) {
			tmp = res;
			res += pre;
			pre = tmp;
		}
		return res;
	}
/****************************test code****************************/
/*
	fibonacciList f(9);

	cout << "getFibonacciNN: " << f.getFibonacciNN() << endl;
	cout << "getFibonacciN: "  << f.getFibonacciN()  << endl;
****debug: getFibonacciNN: 34
           getFibonacciN: 34
*/
/*****************************end test****************************/
};
//#11
class stepPath {
private:
	int N;
	int num;
private:
	int getPathNum(int _n) {
		if (_n < 1)
			return 0;
		if (_n == 1 || _n == 2)
			return _n;
		return getPathNum(_n - 1) + getPathNum(_n - 2);
	}
public:
	stepPath(int _N = 0) : N(_N), num(0) {}
	int getPathNum(void) {
		return getPathNum(N);
	}
};
//#12
class minCoins {
private:
	vector<int> arr;
	int sum;
private:
	int getExchangeNumRec(int _index, int _remain) {
		int res = 0;
		set<int> tmp;

		if (_index == arr.size())
			res = (_remain == 0 ? 1 : 0);
		else{
			for (int i = 0; arr[_index] * i <= _remain; ++i)
				res += getExchangeNumRec(_index + 1, _remain - arr[_index] * i);
		}
		return res;
	}
	int getMinExchangeRes(int _index, int _remain, set<int> &x) {
		int res = 0;
		int tmp = 0;
		int step = 0;

		if (_index == arr.size())
			res = (_remain == 0 ? 1 : 0);
		else {
			for (int i = 0; arr[_index] * i <= _remain; ++i) {
				tmp = getMinExchangeRes(_index + 1, _remain - arr[_index] * i, x);
				if (tmp != 0) {
					x.insert(tmp);
					res += tmp;
				}
			}
		}
		return res;
	}
public:
	minCoins(vector<int> _arr, int _sum = 0) : arr(_arr), sum(_sum) {}
	int getExchangeNumRec(void) {
		return getExchangeNumRec(0, sum);
	}

	int getMinExchangeRec(void) {
		set<int> tmp;
		getMinExchangeRes(0, sum, tmp);
		return *tmp.begin();
	}
	int getMinExchangeDP(void) {
		if (arr.size() == 0 || sum < 0)
			return -1;
		vector< vector<int> > dp(arr.size(), vector<int>(sum + 1));
		for (int i = 1; i <= sum; ++i) {
			dp[0][i] = INT_MAX;
			if ((i - arr[0]) >= 0 && (dp[0][i - arr[0]] != INT_MAX))
				dp[0][i] = dp[0][i - arr[0]] + 1;
		}
		int lchild_up = 0;
		for (size_t i = 1; i < arr.size(); ++i) {
			for (int j = 1; j <= sum; ++j) {
				lchild_up = INT_MAX;
				if (((j - arr[i]) >= 0) && (dp[i][j - arr[i]] != INT_MAX))
					lchild_up = dp[i][j - arr[i]] + 1;
				dp[i][j] = MIN(lchild_up, dp[i - 1][j]);
			}
		}
		return dp[arr.size() - 1][sum] != INT_MAX ? dp[arr.size() - 1][sum] : -1;
	}
	int getMinExchangeDPRes(void) {
		if (sum == 0 || arr.size() == 0)
			return 0;
		int *dp = new int[sum + 1];
		int lchild_up = 0;
		memset(dp, 0, sizeof(int) * (sum + 1));
		for (int i = 1; i <= sum; ++i) {
			dp[i] = INT_MAX;
			if ((i - arr[0] >= 0) && (dp[i - arr[0]] != INT_MAX))
				dp[i] = dp[i - arr[0]] + 1;
		}
		for (size_t i = 1; i < arr.size(); ++i) {
			for (int j = 1; j <= sum; ++j) {
				lchild_up = INT_MAX;
				if ((j - arr[i] >= 0) && (dp[j - arr[i]] != INT_MAX))
					lchild_up = dp[j - arr[i]] + 1;
				dp[j] = MIN(dp[j], lchild_up);
			}
		}
		lchild_up = dp[sum];
		delete[] dp;
		return lchild_up != INT_MAX ? lchild_up : -1;
	}
	int getFixMinExchangeDP(void) {
		if ((arr.size() == 0) || sum < 0)
			return 0;
		vector< vector<int> > dp(arr.size(), vector<int>(sum + 1));
		for (int i = 1; i <= sum; ++i)
			dp[0][i] = INT_MAX;
		if (arr[0] < sum)
			dp[0][arr[0]] = 1;
		int lchild_up = 0;
		for (size_t i = 1; i < arr.size(); ++i) {
			for (int j = 1; j <= sum; ++j) {
				lchild_up = INT_MAX;
				if ((j - arr[i] >= 0) && (dp[i - 1][j - arr[i]] != INT_MAX))
					lchild_up = dp[i - 1][j - arr[i]] + 1;
				dp[i][j] = MIN(dp[i - 1][j], lchild_up);
			}
		}
		return dp[arr.size() - 1][sum] != INT_MAX ? dp[arr.size() - 1][sum] : -1;
	}
	int getFixMinExchangeDPRes(void) {
		if ((arr.size() == 0) || (sum == 0))
			return 0;

		int lchild_up = 0;
		vector<int> dp(sum + 1);
		for (int i = 1; i <= sum; ++i)
			dp[i] = INT_MAX;
		if (arr[0] <= sum)
			dp[arr[0]] = 1;
		for (size_t i = 1; i < arr.size(); ++i) {
			for (int j = sum; j > 0; --j) {
				lchild_up = INT_MAX;
				if ((j - arr[i] >= 0) && (dp[j - arr[i]] != INT_MAX))
					lchild_up = dp[j - arr[i]] + 1;
				dp[j] = MIN(lchild_up, dp[j]);
			}
		}
		return dp[sum] == INT_MAX ? -1 : dp[sum];
	}
/****************************test code****************************/
/*
	vector<int> src;
	src.push_back(1);
	src.push_back(2);
	//src.push_back(6);
	//src.push_back(10);
	minCoins m(src, 3);

	cout << "getMinExchangeRec: " << m.getMinExchangeRec() << endl;
	cout << "getMinExchangeDPRes: " << m.getMinExchangeDPRes() << endl;
	//cout << "getFixMinExchangeDP: " << m.getFixMinExchangeDP() << endl;
	//cout << "getFixMinExchangeDPRes: " << m.getFixMinExchangeDPRes() << endl;
****debug: getFibonacciNN: 34
	getFibonacciN: 34
*/
/*****************************end test****************************/
};
//#13
class KMP {
private:
	string str;
	string match;
private:
	vector<int> getPreFixArr(void) {
		size_t cn = 0, pos = 2;
		vector<int> tmp(match.length(), -1);

		if (match.length() == 1)
			return tmp;
		tmp[0] = -1;
		tmp[1] = 0;
		while (pos < match.length()) {
			if (match[pos - 1] == match[cn])
				tmp[pos++] = tmp[cn++] + 1;    //or tmp[pos++] = ++cn;
			else if (cn > 0)
				cn = tmp[cn];
			else
				tmp[pos++] = 0;
		}
		return tmp;
	}
public:
	KMP(string &_str, string &_match) : str(_str), match(_match) {}
	int getIndexOfByNormal(void) {
		size_t index_s = 0, index_m = 0, tmp = 0;

		while ((index_m < match.length()) && (index_s < str.length())) {
			if (match[index_m] == str[index_s])
				++index_m;
			else {
				index_m = 0;
				index_s = tmp++;
			}
			++index_s;
		}
		if (index_m == match.length())
			return (index_s - match.length());
		return -1;
	}
	int getIndexOfByKMP(void) {
		size_t index_s = 0, index_m = 0, tmp = 0;

		vector<int> prefixarr = getPreFixArr();

		while ((index_s < str.length()) && (index_m < match.length())) {
			if (match[index_m] == str[index_s]) {
				++index_m;
				++index_s;
			}
			else if (prefixarr[index_m] == -1)
				++index_s;
			else
				index_m = prefixarr[index_m];
		}
		return (index_m == match.length() ? index_s - index_m : -1);
	}
/****************************test code****************************/
/*
	string src = "abababc";
	string match = "abc";
	KMP k(src, match);
	cout << "getIndexOfByNormal: " << k.getIndexOfByNormal() << endl;
	cout << "getIndexOfByKMP: " << k.getIndexOfByKMP() << endl;
****debug: 6 9 1 3 1 2 2 5 6 1 1 5 9 7 2 5 6 1 9
	getMinKNumByHeap: 5 3 2 2 1 1 2 1 1 1
	getMinKNumByBFPRT: 1 3 1 2 2 1 1 2 1 5
*/
/*****************************end test****************************/
};
//#14
class Manacher {
private:
	string src;
	string expand;
private:
	void expandString(void) {
		for (size_t i = 0; i < src.length(); ++i) {
			expand += src[i];
			expand += '#';
		}
	}
public:
	Manacher(string &_src) : src(_src), expand("#") {}
	int getMaxLcpsLengthNormal(void) {
		int max_len = 0;

		if (src.length() == 0)
			return 0;
		expandString();
		for (int i = 0; i < expand.length(); ++i) {
			int j = 0;
			while ((i - j >= 0) && (i + j < expand.length())) {
				if (expand[i - j] == expand[i + j])
					++j;
				else 
					break;
			}
			max_len = MAX((2 * j - 1), max_len);
		}
		return (max_len >> 1);
	}
	int getMaxLcpsLengthManacher(void) {
		if (src.length() == 0)
			return 0;

		expand.clear();
		expand += '#';
		expandString();
		vector<int> pArr(expand.length(), 0);
		int index = -1;
		int pR = -1;
		int max_len = INT_MIN;

		for (size_t i = 0; i != expand.length(); ++i) {
			pArr[i] = pR > i ? MIN(pArr[(index << 1) - i], pR - i) : 1;
			while ((i - pArr[i] > -1) && (i + pArr[i] < expand.length())) {
				if (expand[i + pArr[i]] == expand[i - pArr[i]])
					++pArr[i];
				else
					break;
			}
			if (i + pArr[i] > pR) {
				pR = i + pArr[i];
				index = i;
			}
			max_len = MAX(pArr[i], max_len);
		}
		return (max_len - 1);
	}
/****************************test code****************************/
/*
	string src = "111311111";
	Manacher s(src);

	cout << "getMaxLcpsLengthNormal: "   << s.getMaxLcpsLengthNormal() << endl;
	cout << "getMaxLcpsLengthManacher: " << s.getMaxLcpsLengthManacher() << endl;
****debug: getMaxLcpsLengthNormal: 7
           getMaxLcpsLengthManacher: 7
*/
/*****************************end test****************************/
};
//#15
class BFPRT {
private:
	vector<int> src;
	size_t k;
private:
	void heapInsert(vector<int> &_heap, int _value, size_t _index) {
		_heap[_index] = _value;
		while (_index != 0) {
			int parent = ((_index - 1) >> 1);
			if (_heap[parent] < _heap[_index]) {
				swap(_heap[parent], _heap[_index]);
				_index = parent;
			}
			else
				break;
		}
	}
	void heapify(vector<int> &_heap, size_t _index, size_t _heapsize) {
		size_t lchild    = (_index << 1) + 1;
		size_t rchild   = (_index << 1) + 2;
		size_t largest = _index;
		while (lchild < _heapsize) {
			if (_heap[lchild] > _heap[_index])
				largest = lchild;
			if ((rchild < _heapsize) && (_heap[rchild] > _heap[largest]))
				largest = rchild;
			if (largest != _index)
				swap(_heap[largest], _heap[_index]);
			else
				break;
			_index = largest;
			lchild  = (_index << 1) + 1;
			rchild = (_index << 1) + 2;
		}
	}

	int select(vector<int> &_tmp, size_t _begin, size_t _end, size_t _i) {
		if (_begin == _end)
			return _tmp[_begin];
		int pivot = medianOfMedians(_tmp, _begin, _end);
		vector<int> pivot_range = partition(_tmp, _begin, _end, pivot);
		if ((_i >= pivot_range[0]) && (_i <= pivot_range[1]))  // k == _i
			return _tmp[_i];
		else if (_i < pivot_range[0])
			return select(_tmp, _begin, pivot_range[0] - 1, _i);
		else
			return select(_tmp, pivot_range[1] + 1, _end, _i);

		return 0;
	}
	vector<int> partition(vector<int> &_tmp, size_t _begin, size_t _end, int _pivotvalue) {
		vector<int> res(2);
		int small = _begin - 1;
		int cur = _begin;
		int big = _end + 1;
		while (cur != big) {
			if (_tmp[cur] < _pivotvalue)
				swap(_tmp[++small], _tmp[cur++]);
			else if (_tmp[cur] > _pivotvalue)
				swap(_tmp[cur], _tmp[--big]);
			else
				++cur;
		}
		res[0] = small + 1;
		res[1] = big - 1;

		return res;
	}
	int medianOfMedians(vector<int> &_tmp, size_t _begin, size_t _end) {
		int num = _end - _begin + 1;
		int offset = (num % 5 == 0 ? 0 : 1);
		vector<int> res(num / 5 + offset, 0);
		for (size_t i = 0; i < res.size(); ++i) {
			int beginI = _begin + i * 5;
			int endI = beginI + 4;
			res[i] = getMedian(_tmp, beginI, MIN(_end, endI));
		}

		return select(res, 0, res.size() - 1, res.size() / 2);
	}
	int getMedian(vector<int> &_tmp, size_t _begin, size_t _end) {
		insertionSort(_tmp, _begin, _end);
		int sum = _begin + _end;
		int mid = (sum / 2) + (sum % 2);
		return _tmp[mid];
	}
	void insertionSort(vector<int> &_tmp, size_t _begin, size_t _end) {
		for (int i = _begin + 1; i != _end + 1; ++i) {
			for (int j = i; j != _begin; --j) {
				if (_tmp[j - 1] > _tmp[j])
					swap(_tmp[j - 1], _tmp[j]);
				else
					break;
			}
		}
	}
public:
	BFPRT(vector<int> &_src, size_t _k) : src(_src), k(_k) {}
	vector<int> getMinKNumByHeap(void) {
		if ((src.size() < k) || (k < 1))
			return src;
		vector<int> heap(k, 0);
		for (size_t i = 0; i != k; ++i)
			heapInsert(heap, src[i], i);
		for (size_t i = k; i != src.size(); ++i) {
			if (src[i] < heap[0]) {
				heap[0] = src[i];
				heapify(heap, 0, k);
			}
		}
		return heap;
	}
	vector<int> getMinKNumByBFPRT(void) {
		if ((k < 1) || (k > src.size()))
			return src;
		vector<int> tmp = src;
		vector<int> res(k, 0);
		int index = 0;

		int minKth = select(tmp, 0, tmp.size() - 1, k - 1);

		for (int i = 0; i != src.size(); ++i) {
			if (src[i] < minKth)
				res[index++] = src[i];
		}
		for (; index != res.size(); ++index)
			res[index] = minKth;

		return res;
	}
/****************************test code****************************/
/*
	int tmp[] = { 6, 9, 1, 3, 1, 2, 2, 5, 6, 1, 1, 5, 9, 7, 2, 5, 6, 1, 9 };
	vector<int> src(tmp, tmp + 19);
	BFPRT arr(src, 10);
	vectorPrint(src);
	vector<int> res = arr.getMinKNumByHeap();
	cout << "getMinKNumByHeap: ";
	vectorPrint(res);

	vector<int> x = arr.getMinKNumByBFPRT();
	cout << "getMinKNumByBFPRT: ";
	vectorPrint(x);

****debug: 6 9 1 3 1 2 2 5 6 1 1 5 9 7 2 5 6 1 9
           getMinKNumByHeap: 5 3 2 2 1 1 2 1 1 1
           getMinKNumByBFPRT: 1 3 1 2 2 1 1 2 1 5
*/
/*****************************end test****************************/
};
//#16
class subArrs {
private:
	vector<int> src;
	int aim;
	int res;
public:
	subArrs(vector<int> _src, int _aim = 0) : src(_src), aim(_aim), res(0) {}
	int getSubArrs(void) {
		list<size_t> max_list, min_list;
		size_t index_l = 0, index_r = 0;
		while (index_l < src.size()) {
			while (index_r < src.size()) {
				while (!min_list.empty() && (src[min_list.back()] >= src[index_r]))
					min_list.pop_back();
				min_list.push_back(index_r);
				while (!max_list.empty() && (src[max_list.back()] <= src[index_r]))
					max_list.pop_back();
				max_list.push_back(index_r);
				if (src[max_list.front()] - src[min_list.front()] > aim)
					break;
				++index_r;
			}
			if (min_list.front() == index_l)
				min_list.pop_front();
			if (max_list.front() == index_l)
				max_list.pop_front();
			res += (index_r - index_l);
			++index_l;
		}
		return res;
	}
/****************************test code****************************/
/*
	vector<int> src({ 3, 9, 6, 1 });
	subArrs x(src, 4);
	cout << "getSubArrs: " << x.getSubArrs() << endl;

****debug: getSubArrs: 5
*/
/*****************************end test****************************/
};
//#17
class maxTree {
private:
	vector<int> src;
	struct x {
		size_t cur_ss;
		size_t lchild_ss;
		size_t rchild_ss;
	};
	struct nod {
		int value;
		nod *lchild, *rchild;
		nod(int _value = 0, nod *_lchild = nullptr, nod *_rchild = nullptr) :
			value(_value), lchild(_lchild), rchild(_rchild) {}
	};
	nod src_root;
private:
	void addNode(nod *_father, nod *_child) {
		_father->lchild == nullptr ? _father->lchild = _child : _father->rchild = _child;
	}
	void levelTraverse(void) {
		queue<nod *> que;
		nod *tmp = nullptr;

		que.push(&src_root);
		cout << src[que.front()->value] << ' ';
		while (!que.empty()) {
			tmp = que.front();
			if (tmp->lchild == nullptr)
				cout << '#';
			else
				cout << src[tmp->lchild->value];
			cout << ' ';
			if (tmp->rchild == nullptr)
				cout << '#';
			else
				cout << src[tmp->rchild->value];
			que.pop();
			if (tmp->lchild  != nullptr)
				que.push(tmp->lchild);
			if (tmp->rchild != nullptr)
				que.push(tmp->rchild);
		}
	}
public:
	maxTree(vector<int> _src) : src(_src) {}
	void generateMaxTree(void) {
		stack<size_t> m;
		vector<x> res;
		x tmp_res;

		/* get lchild and rchila */
		for (size_t i = 0; i < src.size(); ++i) {
			if (m.empty() == true || src[m.top()] > src[i])
				m.push(i);
			else {
				tmp_res.cur_ss = m.top();
				tmp_res.rchild_ss = i;
				m.pop();
				tmp_res.lchild_ss = (m.empty() ? -1 : m.top());
				res.push_back(tmp_res);
				m.push(i);
			}
		}
		while (!m.empty()) {
			tmp_res.cur_ss = m.top();
			tmp_res.rchild_ss = -1;
			m.pop();
			tmp_res.lchild_ss = (m.empty() ? -1 : m.top());
			res.push_back(tmp_res);
		}

		/* generate binary tree */
		vector<nod> tree;
		for (size_t i = 0; i < src.size(); ++i)
			tree.push_back(*(new nod(i)));
		for (size_t i = 0; i < src.size(); ++i) {
			if (res[i].lchild_ss != -1 && res[i].rchild_ss != -1) {
				if (src[res[i].lchild_ss] < src[res[i].rchild_ss])
					addNode(&tree[res[i].lchild_ss], &(tree[res[i].cur_ss]));
				else
					addNode(&tree[res[i].rchild_ss], &(tree[res[i].cur_ss]));
			}
			else if(res[i].lchild_ss == -1 && res[i].rchild_ss != -1)
				addNode(&tree[res[i].rchild_ss], &(tree[res[i].cur_ss]));
			else if(res[i].lchild_ss != -1 && res[i].rchild_ss == -1)
				addNode(&tree[res[i].lchild_ss], &(tree[res[i].cur_ss]));
		}
		src_root = tree[res.back().cur_ss];

		levelTraverse();  //traverse tree by level
	}
};
//#18 
class maxSubMatrix {
private:
	vector< vector<int> > src;
	int max_s;
public:
	maxSubMatrix(vector< vector<int> > _src) : src(_src), max_s(0) {}
	int getMaxSubMatrixByNormal(void) {
		vector<int> length(src[0].size());
	    int tmp_s;
		int cur_len = 0;
		int index_l = 0;

		for (int col = 0; col < src[0].size(); ++col) {
			for (int row = 0; row < src.size(); ++row) 
					length[col] = (src[row][col] == 1 ? ++length[col] : 0);
		}
		
		for (size_t index_cur = 0; index_cur != length.size(); ++index_cur) {
			if (length[index_cur] != 0) {
				tmp_s = length[index_cur];
				for (int index_l = index_cur - 1; index_l != -1; --index_l) {
					if (length[index_l] >= length[index_cur])
						tmp_s += length[index_cur];
					else
						break;
				}
				for (int index_r = index_cur + 1; index_r != length.size(); ++index_r) {
					if (length[index_r] >= length[index_cur])
						tmp_s += length[index_cur];
					else 
						break;
				}
			}
			max_s = MAX(max_s, tmp_s);
		}
		return max_s;
	}
	int getMaxSubMatrixByStack(void) {
		vector<int> length(src[0].size());
		stack<int> s;
		int tmp_s = 0;

		for (size_t col = 0; col < src[0].size(); ++col) {
			for (size_t row = 0; row < src.size(); ++row)
				length[col] = (src[row][col] == 1 ? ++length[col] : 0);
		}
		s.push(0);
		for (size_t index_cur = 1; index_cur != src[0].size(); ++index_cur) {
			if (length[index_cur] != 0) {
				if (length[index_cur] >= length[s.top()])
					s.push(index_cur);
				else {
					while (!s.empty() && length[s.top()] > length[index_cur]) {
						tmp_s = (index_cur - s.top()) * length[s.top()];
						max_s = MAX(max_s, tmp_s);
						s.pop();
					}
					s.push(index_cur);
				}
			}
		}
		while (!s.empty()) {
			tmp_s = (src[0].size() - (s.empty() == true ? 0 : s.top())) * length[s.top()];
			max_s = MAX(max_s, tmp_s);
			s.pop();
		}

		return max_s;
	}
/****************************test code****************************/
/*
	vector< vector<int> > src = { {1, 0, 1, 1},
	                              {1, 1, 1, 1},
								  {1, 1, 1, 0},
								  {1, 1, 1, 1} };
	maxSubMatrix x(src);
	cout << "getMaxSubMatrixByNormal: " << x.getMaxSubMatrixByNormal() << endl;
	cout << "getMaxSubMatrixByStack: "  << x.getMaxSubMatrixByStack() << endl;
****debug: getMaxSubMatrixByNormal: 9
           getMaxSubMatrixByStack: 9
*/
/*****************************end test****************************/
};
//#19
class securityPlan {
private:
	vector<int> src;
	int max_path;
	struct element {
		int times;
		int script;
		element(int _script = 0, int _times = 1) : times(_times), script(_script) {}
	};
private:
	void push_stack(stack<element> &_s, int _cur) {
		if (src[_cur] == src[_s.top().script])
			++_s.top().times;
		else
			_s.push(_cur);
	}
public:
	securityPlan(vector<int> _src) : src(_src), max_path(0) {}
	int getSecurityPlan(void) {
#define GETNEXTSCRIPT(cur)  ((cur + 1) < src.size() ? (cur + 1) : 0)
#define GETCOMBINATION(x)   (x == 1 ? 0 : ((x * (x - 1)) >> 1))
		stack<element> s;
		element tmp;
		size_t max_s = 0;
		size_t cur_s = 0;

		for (size_t i = 1; i < src.size(); ++i)
			max_s = src[max_s] > src[i] ? max_s : i;
		s.push(max_s);
		cur_s = GETNEXTSCRIPT(max_s);
		for (size_t times = 1; times < src.size(); ++times) {
			if (src[s.top().script] >= src[cur_s])
				push_stack(s, cur_s);
			else {
				while (src[s.top().script] < src[cur_s]) {
					tmp = s.top();
					s.pop();
					max_path += GETCOMBINATION(tmp.times) + (tmp.times << 1);
				}
				if (src[s.top().script] >= src[cur_s])
					push_stack(s, cur_s);
			}

			cur_s = GETNEXTSCRIPT(cur_s);
		}
		while (s.size() > 1) {
			tmp = s.top();
			s.pop();
			max_path += GETCOMBINATION(tmp.times) + (s.top().times > 1 ? (tmp.times << 1) : tmp.times);
		}
		max_path += GETCOMBINATION(s.top().times);

		return max_path;
	}
	/****************************test code****************************/
	/*
	vector<int> src({3, 3, 3, 2, 2, 3, 3, 3, 2, 3});
	securityPlan x(src);
	cout << "getSecurityPlan: " << x.getSecurityPlan() << endl;
	****debug: getSecurityPlan: 28
	*/
	/*****************************end test****************************/
};
//#20
class morris {
private:
	struct nod {
		int value;
		nod *lchild;
		nod *rchild;
		nod(int _value = 0, nod *_lchild = nullptr , nod *_rchild = nullptr) :
			value(_value), lchild(_lchild), rchild(_rchild) {}
	};
	nod *src_root;
private:
	void printRightPath(nod *_head) {
		stack<int> tmp;
		nod *most_rchild = _head->lchild;

		while (most_rchild != nullptr && most_rchild != _head) {
			tmp.push(most_rchild->value);
			most_rchild = most_rchild->rchild;
		}
		stackPrint(tmp);
	}
	void printRightBoundary(void) {
		stack<int> x;
		nod *tmp = src_root;

		while (tmp != nullptr) {
			x.push(tmp->value);
			tmp = tmp->rchild;
		}
		stackPrint(x);
	}
public:
	morris(nod *_root = nullptr) : src_root(_root) {}
	void generateTree(void) {
		queue<nod *> que;
		nod *tmp = nullptr;
		int ldata = 0, rdata = 0, x = 0;

		cout << "Please enter root node's value: ";
		cin >> x;
		src_root = new nod(x);
		que.push(src_root);
		while (!que.empty()) {
			tmp = que.front();
			que.pop();
			cout << "plaese enter " << tmp->value << "'s ldata and rdata: ";
			cin >> ldata >> rdata;
			if (ldata != 0)
				que.push(tmp->lchild = new nod(ldata));
			if (rdata != 0)
				que.push(tmp->rchild = new nod(rdata));
		}
	}
	void morrisPre(void) {
		nod *cur = nullptr, *most_rchild = nullptr;
	    
		cur = src_root;
		while (cur != nullptr) {
			most_rchild = cur->lchild;
			if (most_rchild != nullptr) {
				while (most_rchild->rchild != nullptr && most_rchild->rchild != cur)
					most_rchild = most_rchild->rchild;
				if (most_rchild->rchild == nullptr) {
					cout << cur->value << ' ';
					most_rchild->rchild = cur;
					cur = cur->lchild;
					continue;
				}
				else
					most_rchild->rchild = nullptr;
			}
			else
				cout << cur->value << ' ';
			cur = cur->rchild;
		}
	
	}
	void morrisMid(void) {
		nod *cur = nullptr, *most_rchild = nullptr;

		cur = src_root;
		while (cur != nullptr) {
			most_rchild = cur->lchild;
			if (most_rchild != nullptr) {
				while (most_rchild->rchild != nullptr && most_rchild->rchild != cur)
					most_rchild = most_rchild->rchild;
				if (most_rchild->rchild == nullptr) {
					most_rchild->rchild = cur;
					cur = cur->lchild;
					continue;
				}
				else
					most_rchild->rchild = nullptr;
			}
			cout << cur->value << ' ';
			cur = cur->rchild;
		}
	}
	void morrisPost(void) {
		nod *cur = nullptr, *most_rchild = nullptr;
		vector<nod *> tmp;

		cur = src_root;
		while (cur != nullptr) {
			most_rchild = cur->lchild;
			if (most_rchild != nullptr) {
				while ((most_rchild->rchild != nullptr) && (most_rchild->rchild != cur))
					most_rchild = most_rchild->rchild;
				if (most_rchild->rchild == nullptr) {
					most_rchild->rchild = cur;
					cur = cur->lchild;
					continue;
				}
				else {
					printRightPath(cur);
					most_rchild->rchild = nullptr;
				}
			}
			cur = cur->rchild;
		}
		printRightBoundary();
	}
/****************************test code****************************/
/*
	morris x;;
	x.generateTree();
	cout << "morrisPre: ";
	x.morrisPre();
	cout << endl;
	cout << "morrisMid: ";
	x.morrisMid();
	cout << endl;
	cout << "morrisPost: ";
	x.morrisPost();
	cout << endl;
****debug: Please enter root node's value: 1
		   plaese enter 1's ldata and rdata: 2 3
		   plaese enter 2's ldata and rdata: 4 5
		   plaese enter 3's ldata and rdata: 6 7
		   plaese enter 4's ldata and rdata: 0 0
		   plaese enter 5's ldata and rdata: 0 0
		   plaese enter 6's ldata and rdata: 0 0
		   plaese enter 7's ldata and rdata: 0 0
		   morrisPre: 1 2 4 5 3 6 7
		   morrisMid: 4 2 5 1 6 3 7
		   morrisPost: 4 5 2 6 7 3 1
*/
/*****************************end test****************************/
};
//#21  continuing...
class skylineProblem {  
private:
	struct nod {
		size_t position;
		int len;
		bool isup;
		nod(size_t _position = 0, int _len = 0, bool _isup = false)
			: position(_position), len(_len), isup(_isup) {}
	};
	vector< vector<int> > src;
	vector< vector<int> > sktline;

	vector< nod *> tmp;
private:

public:
	skylineProblem(vector< vector<int> > _src) : src(_src) {}
	vector< vector<int> > getSkyline(void) {
	    

	}
};  
//#22
class maxSubArrLenAsSum {
private:
	vector<int> src;
	int aim;
	int len;
public:
	maxSubArrLenAsSum(vector<int> _src, int _aim) : src(_src), aim(_aim), len(-1) {}
	int getMaxSubArrLen1(void) {
		int tmp = src[0];
		int l_script = 0, r_script = 0;

		while(r_script < src.size()) {
			if (tmp < aim) {
				++r_script;
				if (r_script == src.size())
					break;
				tmp += src[r_script];
			}
			else if (tmp == aim) {
				len = MAX(len, (r_script - l_script + 1));
				tmp -= src[l_script++];
			}
			else
				tmp -= src[l_script++];
		}

		return len;
	}
	int getMaxSubArrLen2(void) {
		unordered_map<int, int> firstPos;
		int sum = 0;

		firstPos.insert(pair<int, int>(0, -1));
		for (int index = 0; index != src.size(); ++index) {
			sum += src[index];
			if (firstPos.count(sum - aim) != 0) 
				len = MAX(index - firstPos.at(sum - aim), len);
			if (firstPos.count(sum) == 0)
				firstPos.insert(pair<int, int>(sum, index));
		}

		return len;
	}
/****************************test code****************************/
/*
	maxSubArrLenAsSum x({ 7, 3, 2, 1, 1, -1, -6, 7, 8 }, 7);
	maxSubArrLenAsSum y({ 7, 3, 2, 1, 1, 7, 8 }, 7);

	cout << "getMaxSubArrLen2: " << x.getMaxSubArrLen2() << endl;
	cout << "getMaxSubArrLen1: " << y.getMaxSubArrLen1() << endl;
****debug: getMaxSubArrLen2: 7
           getMaxSubArrLen2: 4
*/
/*****************************end test****************************/
};
//#23
class mostxSubXORArrNum {
private:
	vector<int> src;
	int num;
public:
	mostxSubXORArrNum(vector<int> _src, int _num = 0) : src(_src), num(_num) {}
	int getMostSubXORArrNum(void) {
		vector<int> mosts(src.size(), 0);
		int xor = 0;
		unordered_map<int, int> tmp;

		tmp.insert(pair<int, int>(0, -1));
		for (size_t index = 0; index < src.size(); ++index) {
			xor ^= src[index];
			if (tmp.count(xor)) {
				int pre = tmp.find(xor)->second;
				mosts[index] = (pre == -1 ? 1 : (mosts[pre] + 1));
			}
			if (index > 0)
				mosts[index] = MAX(mosts[index - 1], mosts[index]);
			if (tmp.count(xor))
				tmp.erase(xor);    //refresh unorder_map
			tmp.insert(pair<int, int>(xor, index));
			num = MAX(num, mosts[index]);
		}

		return num;
	}
	int comparator(void) {
		if (src.size() == 0)
			return 0;

		vector<int> xors(src.size());
		int xor = 0;
		for (size_t index = 0; index < src.size(); ++index) {
			xor ^= src[index];
			xors[index] = xor;
		}

		vector<int> mosts(src.size());
		mosts[0] = (src[0] == 0 ? 1 : 0);
		for (size_t i = 1; i < src.size(); ++i) {
			mosts[i] = (xors[i] == 0 ? 1 : 0);
			for (size_t j = 0; j < i; ++j) {
				if ((xors[i] ^ xors[j]) == 0)
					mosts[i] = MAX(mosts[i], mosts[j] + 1);
			}
			mosts[i] = MAX(mosts[i], mosts[i - 1]);
		}
	
		return mosts[mosts.size() - 1];
	}
/****************************test code****************************/
/*
	vector<int> arr;
	mostxSubXORArrNum x({ 1, 2, 3, 0, 1, 2, 3 });
	cout << "getMostSubXORArrNum: " << x.getMostSubXORArrNum() << endl;
	cout << "comparator: " << x.comparator() << endl;

	for (int time = 0; time < testTime; ++time) {
		//cout << time << endl;
		arr = generateRandomArray(maxSize, maxValue); 
		mostxSubXORArrNum y(arr);
		int res = y.getMostSubXORArrNum();
		int tmp = y.comparator();
		if (res != tmp) {
			succeed = false;
			vectorPrint<int>(arr);
			cout << "times: " << time << endl;
			cout << "getMostSubXORArrNum: " << res << endl;
			cout << "comparator: " << tmp << endl;
			break;
		}
	}
****debug: getMostSubXORArrNum: 3
           comparator: 3
           Nice!
*/
/*****************************end test****************************/
};
//#24
class maxSearchSubTree {
private:

	nod *src_root;
	struct returnType {
		int tree_size;
		nod *head;
		int max_value;
		int min_value;
		returnType(int _tree_size = 0, nod *_head = nullptr, int _max_value = INT_MIN, int _min_value = INT_MAX) :
			tree_size(_tree_size), head(_head), max_value(_max_value), min_value(_min_value) {}
	};
	returnType res;
private:
	returnType getMaxSearchSubTreeRec(nod *_head) {
		returnType res(0, nullptr, INT_MIN, INT_MAX);

		if (_head == nullptr) 
			return res;
		else {
			returnType lchild;
			returnType rchild;

			lchild = getMaxSearchSubTreeRec(_head->lchild);
			rchild = getMaxSearchSubTreeRec(_head->rchild);
			//#0 lchild + _head + rchild
			if (lchild.head == _head->lchild &&
				rchild.head == _head->rchild &&
				_head->value > lchild.max_value &&
				_head->value < rchild.min_value) {
				res.tree_size = lchild.tree_size + 1 + rchild.tree_size;
				res.head = _head;
				res.max_value = MAX(_head->value, lchild.max_value);
				res.min_value = MIN(_head->value, rchild.min_value);
				return res;
			}
			else  //max(lchild.size, rchild.size)
				return (lchild.tree_size > rchild.tree_size ? lchild : rchild);
		}
	}
	returnType comparator(nod *_head) {
		return res;
		if (_head == nullptr) {
			return res;
		}
		//node *lchild = _head->lchild;
		returnType lchildSubTressInfo = comparator(_head->lchild);
		//node *rchild = _head->rchild;
		returnType rchildSubTressInfo = comparator(_head->rchild);

		int includeItSelf = 0;
		if (lchildSubTressInfo.head == _head->lchild
			&&rchildSubTressInfo.head == _head->rchild
			&& _head->value > lchildSubTressInfo.max_value
			&& _head->value < rchildSubTressInfo.min_value
			) {
			includeItSelf = lchildSubTressInfo.tree_size + 1 + rchildSubTressInfo.tree_size;
		}
		int p1 = lchildSubTressInfo.tree_size;
		int p2 = rchildSubTressInfo.tree_size;
		int maxSize = MAX(MAX(p1, p2), includeItSelf);

		nod *maxHead = (p1 > p2 ? lchildSubTressInfo.head : rchildSubTressInfo.head);
		if (maxSize == includeItSelf) {
			maxHead = _head;
		}

		res.head = maxHead;
		res.tree_size = maxSize;
		res.min_value = MIN(MAX(lchildSubTressInfo.min_value, rchildSubTressInfo.min_value), _head->value);
		res.max_value = MAX(MAX(lchildSubTressInfo.max_value, rchildSubTressInfo.max_value), _head->value);
		return  res;
	}
public:
	maxSearchSubTree(nod *_root = nullptr) : src_root(_root), res(0, nullptr, INT_MAX, INT_MIN) {}
	void generateTree(void) {
		queue<nod *>src;
		nod *tmp = nullptr;
		int lvalue = 0, rvalue = 0, x = 0;

		cout << "please enter root node value: ";
		cin >> x;
		src_root = new nod(x);
		src.push(src_root);
		while (!src.empty()) {
			tmp = src.front();
			src.pop();
			cout << "please enter " << tmp->value << "'s lvalue and rvalue: ";
			cin >> lvalue >> rvalue;
			if (lvalue != -1)
				src.push(tmp->lchild = new nod(lvalue));
			if (rvalue != -1)
				src.push(tmp->rchild = new nod(rvalue));
		}
	}
	int getMaxSearchSubTreeRec(void) {
		res = getMaxSearchSubTreeRec(src_root);

		return res.tree_size;
	}
	int comparator(void) {
		returnType res;
		res = comparator();

		return res.tree_size;
	}
/****************************test code****************************/
/*
	int res = 0;
	int max = 0;
	for (int time = 0; time < test_ime; ++time) {
		node *root = generateRandomBinaryTree(max_size, max_value);
		BinaryTreePrint(root);
		cout << endl;
		maxSearchSubTree x(root);
		res = x.getMaxSearchSubTree();
		max = MAX(max, res);
		cout << "getMaxSearchSubTree: " << res << endl;
	}
	cout << max << endl;
****debug:
*/
/*****************************end test****************************/
};
//#25
class LRU {
private:
	struct nod {
		string key;
		int value;
		nod(string _key, int _value) : key(_key), value(_value) {}
	};
	unordered_map<string, nod*> hash_map;
	list<nod *> dual_list;
	int catch_size;
private:
	void addNode2End(nod *_node) {
		dual_list.push_back(_node);
	}
	void moveNode2End(nod *_node) {
		nod *tmp = dual_list.front();
		dual_list.pop_front();
		dual_list.push_back(tmp);
	}
public:
	LRU(int _catch_size) : catch_size(_catch_size) {}
	void set(string _key, int _value) {
		nod *_node = new nod(_key, _value);
		if (get(_node->key) == INT_MIN && hash_map.size() < catch_size)
		    addNode2End(_node);
		else if (get(_node->key) != INT_MIN){
			hash_map[_node->key]->value = _node->value;
			moveNode2End(_node);
		}
	}
	int get(string _key) {
		if (hash_map.count(_key) != 0) {
			moveNode2End(hash_map.find(_key)->second);
			return hash_map[_key]->value;
		}	
		else
			return INT_MIN;
	}
	void mostFrequent(void) {
		cout << dual_list.back()->key << "->";
cout << dual_list.back()->value << endl;
	}
	/****************************test code****************************/
	/*
		LRU x(3);
		x.set("A", 1);
		x.set("B", 4);
		x.set("C", 2);
		x.set("D", 3);
		x.mostFrequent();
		x.set("B", 5);
		x.mostFrequent();

		cout << (succeed ? "Nice!" : "Fucking fucked!") << endl;
	****debug:D->3
			  B->5
			  Nice!
	*/
	/*****************************end test****************************/
};
//#26  
class LFU {
private:
	struct node_times {
		int times;
		list<string> *tasks;
		node_times(int _times = 0, list<string> *_tasks = nullptr) : times(_times), tasks(_tasks) {}
	};
	list< node_times > times;
	unordered_map<string, int> task_info;  //task's name and operate times
	unordered_map<int, node_times > times_info;  //operate times and tasks
private:
	void addTimesNode(int _times) {
		node_times x(_times, new list<string>);
		times_info.insert(pair<int, node_times>(_times, x));
	}
public:
	LFU() {
		addTimesNode(1);
	}
	void set(string _task) {
		if (task_info.count(_task) == 0) {
			task_info.insert(pair<string, int>(_task, 1));
			times_info[1].tasks->push_back(_task);
		}
		else {
			int times = task_info[_task];
			++task_info[_task];
			node_times tmp = times_info[times];
			tmp.tasks->remove(_task);
			if (times_info.count(times + 1) == 0) {
				addTimesNode(times + 1);
				times_info[times + 1].tasks->push_back(_task);
			}
			else
				times_info[times].tasks->push_back(_task);
		}
		cout << "set/get " << _task << endl;
	}
	int get(string _task) {
		int times = task_info.count(_task);
		set(_task);
		return times + 1;
	}
	void mostOperate(void) {
		string task = (--times_info.end())->second.tasks->back();
		cout << "mostOperateTask: " << task << "    times: " << task_info[task] << endl;
	}
	/****************************test code****************************/
	/*
		LFU x;
		x.set("A");
		x.set("B");
		x.set("C");
		x.get("B");
		x.mostOperate();
	****debug:  set/get A
				set/get B
				set/get C
				set/get B
				mostOperateTask: B    times: 2
	*/
	/*****************************end test****************************/
};
//#27
class bstTopoSBT {
private:
	vector<int> src_node;
	vector<int> res_node;
	nod *src_root;
	nod *res_root;
	struct returnData {
		nod *head;
		int maxValue;
		int minValue;
		returnData(nod *_head = nullptr, int _maxValue = INT_MIN, int _minValue = INT_MAX) :
			head(_head), maxValue(_maxValue), minValue(_minValue) {}
	};
public:
	bstTopoSBT(vector<int> _src_node) : src_node(_src_node), src_root(nullptr), res_root(nullptr) {}
	bstTopoSBT(nod *_src_root = nullptr) : src_root(_src_root), res_root(nullptr) {}
	returnData getBstTopoSBT(nod *_head) {
		returnData  tmp, l_data, r_data;

		if (_head->lchild == nullptr && _head->rchild == nullptr) {
			tmp.head = _head;
			tmp.maxValue = _head->value;
			tmp.minValue = _head->value;
			return tmp;
		}

		l_data = getBstTopoSBT(_head->lchild);
		r_data = getBstTopoSBT(_head->rchild);

		_head->lchild = l_data.head;
		tmp.minValue = l_data.maxValue;
		_head->rchild = r_data.head;
		tmp.maxValue = r_data.maxValue;
		tmp.head = _head;

		if (l_data.maxValue > _head->value) {
			_head->lchild->lchild = nullptr;
			_head->lchild->rchild = nullptr;
			tmp.minValue = _head->lchild->value;
		}
		if (r_data.maxValue < _head->value) {
			_head->rchild->lchild = nullptr;
			_head->rchild->rchild = nullptr;
			tmp.maxValue = _head->rchild->value;
		}
		if (_head->lchild->value > _head->value)
			_head->lchild = nullptr;
		if(_head->rchild->value < _head->value)
			_head->rchild = nullptr;

		return tmp;
	}
	nod *getBstTopoSBT(void) {
		returnData tmp = getBstTopoSBT(src_root);
		return tmp.head;
	}









	/****************************test code****************************/
	/*

	****debug:
	*/
	/*****************************end test****************************/
};
//#28
class binaryTree {
private:
	struct returnData {
		int l_pos;
		int r_pos;
		returnData(int _l_most_r = INT_MAX, int _r_most_l = INT_MIN) : l_pos(_l_most_r), r_pos(_r_most_l) {}
	};
	nod *src_root;
	nod *res_root;
private:
	void nodesPosAdj(nod *_head, int _changes) {
		if (_head == nullptr)
			return;
		if (_head->lchild)
			_head->lchild->pos += _changes;
		if (_head->rchild)
			_head->rchild->pos += _changes;
		nodesPosAdj(_head->lchild, _changes);
		nodesPosAdj(_head->rchild, _changes);
	}

	returnData binaryPosAdj(nod *_root, nod *_head) {
		int changes = 0;
		returnData tmp, lchild_data, rchild_data;

		if (_head == nullptr)
			return tmp;
		lchild_data = binaryPosAdj(_root, _head->lchild);
		rchild_data = binaryPosAdj(_root, _head->rchild);
		while (_head->lchild != nullptr && lchild_data.r_pos >= _head->pos) {
			changes = lchild_data.r_pos - _head->pos + 2;
			_head->pos += changes;
			if (_head->rchild) {
				_head->rchild->pos += changes;
				nodesPosAdj(_head->rchild, changes);
				if (_head != _root)
					rchild_data = binaryPosAdj(_root, _head->rchild);
			}
			else
				rchild_data.l_pos = MIN(rchild_data.l_pos, _head->pos);
		}
		while (_head->rchild != nullptr && rchild_data.l_pos <= _head->pos) {
			changes = _head->pos - rchild_data.l_pos + 2;
			_head->rchild->pos += changes;
			nodesPosAdj(_head->rchild, changes);
			//if (_head != _root)
			rchild_data = binaryPosAdj(_root, _head->rchild);
		}
		tmp.l_pos = MIN(_head->pos, lchild_data.l_pos);
		tmp.r_pos = MAX(_head->pos, rchild_data.r_pos);
		return tmp;
	}

	void refreshNodesPos(int _min_pos, nod *_head) {
		if (_head == nullptr)
			return;
		_head->pos += _min_pos;
		if (_head->lchild)
			refreshNodesPos(_min_pos, _head->lchild);
		if (_head->rchild)
			refreshNodesPos(_min_pos, _head->rchild);
	}

	nod *preTraverse(nod *_head) {
		if (_head == nullptr)
			return nullptr;
		cout << _head->value << ' ';
		preTraverse(_head->lchild);
		preTraverse(_head->rchild);
	}

	nod *midTraverse(nod *_head) {
		if (_head == nullptr)
			return nullptr;
		midTraverse(_head->lchild);
		cout << _head->value << ' ';
		midTraverse(_head->rchild);
	}

	nod *postTraverse(nod *_head) {
		if (_head == nullptr)
			return nullptr;
		postTraverse(_head->lchild);
		postTraverse(_head->rchild);
		cout << _head->value << ' ';
	}
public:
	binaryTree() : src_root(nullptr), res_root(nullptr) {}

	nod* generateFixedBinaryTree(void) {
		queue<nod *> src;
		nod *tmp = nullptr;
		int x = 0, lchild = 0, rchild = 0;

		cout << "please enter root value: ";
		cin >> x;
		src.push(src_root = new nod(x));
		while (!src.empty()) {
			tmp = src.front();
			src.pop();
			cout << "please enter " << tmp->value << "'s lchild and rchild: ";
			cin >> lchild >> rchild;
			if (lchild != -1)
				src.push(tmp->lchild = new nod(lchild));
			if (rchild != -1)
				src.push(tmp->rchild = new nod(rchild));
		}
		return src_root;
	}

	nod *generateFixedBinaryTree(vector<int> _src_nodes) {
		nod *tmp = nullptr;
		queue<nod *> x;
		int index = 0;
		int min_pos = 0;

		if (_src_nodes[0] == '#')
			return nullptr;
		x.push(src_root = new nod(_src_nodes[index++], 0));
		min_pos = src_root->pos;
		while (!x.empty()) {
			tmp = x.front();
			x.pop();
			if (_src_nodes[index] != '#')
				x.push(tmp->lchild = new nod(_src_nodes[index], (tmp->pos - intlen(_src_nodes[index]) - 1)));
			++index;
			if (_src_nodes[index] != '#')
				x.push(tmp->rchild = new nod(_src_nodes[index], (tmp->pos + (intlen(tmp->value) + 1))));
			++index;
		}
		return src_root;
	}

	nod *generateRandomBinaryTree(int _maxsize, int _maxvalue) {
		nod *tmp = nullptr;
		queue< nod * > src;
		int nullvalue = (_maxvalue >> 1);
		int value = 0;
		int maxsize = (rand() % _maxsize);
		int size = 0;

		src_root = new nod(rand() % _maxvalue);
		if (src_root->value == nullvalue)
			return nullptr;
		++size;
		src.push(src_root);
		while ((!src.empty()) && (size < maxsize)) {
			tmp = src.front();
			src.pop();
			if (tmp->value != '#') {
				value = rand() % _maxvalue;
				if (value != nullvalue) {//((value > (nullvalue - 10)) && (value > (nullvalue + 10))) {
					++size;
					src.push(tmp->lchild = new nod(value));
				}
				else
					src.push(tmp->lchild = new nod('#'));
				value = rand() % _maxvalue;
				if (value != nullvalue) {//((value > (nullvalue - 10)) && (value > (nullvalue + 10))) {
					++size;
					src.push(tmp->rchild = new nod(value));
				}
				else
					src.push(tmp->rchild = new nod('#'));
			}
		}
		return src_root;
	}

	void binaryTreePrint(void) {
		queue<nod *>src;
		nod *tmp = nullptr;
		int level = 0;
		int pre_num_pos = 0, pre_line_pos = 0, min_pos = 0, pre_len = 0, fill_len = 0;
		vector<int> nodes_per_level(2, 0);  //level nodes;
		vector<string> chars(1);
		vector<string> nums(1);

		if (src_root == nullptr)
			return;
		min_pos = binaryPosAdj(src_root, src_root).l_pos;
		refreshNodesPos(abs(min_pos), src_root);
		src.push(src_root);
		nodes_per_level[level] = 1;
		while (!src.empty()) {
			tmp = src.front();
			src.pop();
			fill_len = tmp->pos - pre_num_pos - pre_len;
			chars[level] += *new string((tmp->pos - pre_line_pos - pre_len), ' ');
			chars[level] += pre_len > 1 ? *new string(pre_len - 1, ' ') : "";
			chars[level] += '|';
			if (tmp->lchild) {
				nums[level] += *new string(tmp->lchild->pos - (tmp->pos - fill_len), ' ');
				nums[level] += *new string(tmp->pos - tmp->lchild->pos, '-');
			}
			else
				nums[level] += *new string(fill_len, ' ');
			nums[level] += int2str(tmp->value);
			pre_num_pos = tmp->pos;
			pre_line_pos = tmp->pos;
			pre_len = intlen(tmp->value);
			if (tmp->rchild) {
				int len = intlen(tmp->rchild->value);
				len = tmp->rchild->pos - tmp->pos - (len > 1 ? len - 1 : 0);
				nums[level] += *new string(len, '-');
				pre_num_pos += len;
			}
			if (tmp->lchild) {
				src.push(tmp->lchild);
				nodes_per_level[level + 1]++;
			}
			if (tmp->rchild) {
				src.push(tmp->rchild);
				nodes_per_level[level + 1]++;
			}
			if (--nodes_per_level[level] == 0) {
				nums[level] += '\n';
				nums.push_back("");
				chars[level] += '\n';
				chars.push_back("");
				nodes_per_level.push_back(0);
				++level;
				pre_num_pos = 0;
				pre_line_pos = 0;
				pre_len = 0;
			}
		}

		int index_nums = 0, index_chars = 0;
		for (size_t index = 0; index < nums.size() - 1; ++index) {
			cout << nums[index_nums++];
			cout << chars[++index_chars];
		}
	}

	void preTraverse(void) {
		preTraverse(src_root);
	}

	void midTaaverse(void) {
		midTraverse(src_root);
	}

	void postTraverse(void) {
		postTraverse(src_root);
	}
/****************************test code****************************/
/*
	vector<int> src_nodes({        6, 
		                    1,                        12, 
		               0,        3,          10,                13, 
		           '#', '#', '#', '#',    4,     14,        20,     16, 
		                                2,  5, 11, 15, '#', '#', '#', '#',
		                            '#', '#', '#', '#', '#', '#', '#', '#' });
	vector<int> src_nodes2({     2, 
		               3,                   4, 
		           5,       6,        3,                   5, 
		       '#', '#', '#', '#',  1,       2,       3,          4, 
		                        '#', '#','#', '#', '#', '#','#', '#' });
	binaryTree x;
	x.generateFixedBinaryTree(src_nodes);
	x.binaryTreePrint();
****debug:  ----6-------------------
            |                       |
		  --1--           ----------12----
		  |   |           |              |
		  0   3       ----10----      ---13--
					  |        |      |     |
				    --4--   ---14--   20    16
				    |   |   |     |
				    2   5   11    15
nice!
*/
/*****************************end test****************************/
};
//#29
class printTreeEdge {
private:
	nod *root = nullptr;
	vector< vector<nod*> > res;
public:
	printTreeEdge(nod *_root = nullptr) : root(_root) {}
	void getTreeEdge1(void) {
		nod *tmp = nullptr;
		queue<nod *> src;
		int level = 0;
		vector<int> level_nodes(2, 0);
		vector<int> last_level;
		queue<nod*> front;
		stack<nod*> back;

		src.push(root);
		level_nodes[level] = 1;
		front.push(root);
		while (!src.empty()) {
			tmp = src.front();
			last_level.push_back(tmp->value);
			src.pop();
			if (tmp->lchild) {
				src.push(tmp->lchild);
				level_nodes[level + 1]++;
			}
			if (tmp->rchild) {
				src.push(tmp->rchild);
				level_nodes[level + 1]++;
			}
			if (--level_nodes[level] == 0) {
				++level;
				level_nodes.push_back(0);
				back.push(tmp);
				if (!src.empty()) {
					front.push(src.front());
					last_level.clear();
				}
			}
		}
		while (!front.empty()) {
			cout << front.front()->value << ' ';
			front.pop();
		}
		for (size_t index = 1; index < last_level.size() - 1; ++index)
			cout << last_level[index] << ' ';
		while (!back.empty()) {
			cout << back.top()->value << ' ';
			back.pop();
			if (back.size() == 1)
				return;
		}
	}
/****************************test code****************************/
/*
    binaryTree tree;
	vector<int> src_nodes({          1, 
		                   2,                                   3, 
		               '#',    4,                    5,               6,
		                  7,        8,          9,        10,      '#', '#', 
		               '#', '#', '#', 11,     12, '#', '#', '#',
		                           13, 14,  15, 16,
		'#', '#', '#', '#', '#', '#', '#', '#' });
	node *root = tree.generateFixedBinaryTree(src_nodes);
	tree.binaryTreePrint();
	printTreeBoard x(root);
	x.getTreeBoard();
	cout << endl << endl << endl;
	node *_root = tree.generateRandomBinaryTree(30, 30);
	tree.binaryTreePrint();
	printTreeBoard y(_root);
	y.getTreeBoard();
****debug:  ----------------1------------------------------
			|                                             |
			2----                                     ----3--
				|                                     |     |
			  --4--                                 --5-    6
			  |   |                                 |   |
			  7   8----                        -----9   10
					   |                       |
					---11--                 ---12--
					|     |                 |     |
					13    14                15    16
			1 2 4 7 11 13 14 15 16 12 10 6 3


				  --------17-----------------
				  |                         |
			  ----4----                   --10-
			  |       |                   |   |
			--29-   --4--                 18  18
			|   |   |   |
			22  14  5   5
			17 4 29 22 14 5 5 18 10 nice!
*/
/*****************************end test****************************/
};
//#30  continue...
class longestPathInTree {
private:
	nod *root;
	nod *res;
	struct returnData {
		int sum;

	};
public:
	longestPathInTree(nod *_root = nullptr) : root(_root), res(nullptr) {}
	void getLongestPath(void) {

	}


	/****************************test code****************************/
	/*

	****debug:
	*/
	/*****************************end test****************************/
};
//#31
class zigzagPrint {
private:
	nod *root;
	vector< vector<int> > src;
private:
	void intuitionPrint(vector< vector<int> > _tmp, int _index) {
		for (size_t index = 0; index < _tmp[_index].size(); ++index)
			cout << _tmp[_index][index] << ' ';
	}
	void reversePrint(vector< vector<int> > _tmp, int _index) {
		for (int index = _tmp[_index].size() - 1; index >= 0; --index)
			cout << _tmp[_index][index] << ' ';
	}
public:
	zigzagPrint(nod *_root = nullptr) : root(_root) {}
	void print(void) {
		nod *tmp = nullptr;
		queue<nod *> que;
		int level = 0;
		vector< int > node_per_level(2);
		
		que.push(root);
		src.push_back({});
		node_per_level[level]++;
		while (!que.empty()) {
			tmp = que.front();
			que.pop();
			src[level].push_back(tmp->value);
			if (tmp->lchild) {
				que.push(tmp->lchild);
				node_per_level[level + 1]++;
			}
			if (tmp->rchild) {
				que.push(tmp->rchild);
				node_per_level[level + 1]++;
			}
			if (--node_per_level[level] == 0) {
				level++;
				node_per_level.push_back(0);
				src.push_back({});
			}
		}
		for (size_t index = 0; index < src.size() - 1; ++index) {
			if (index % 2 == 0) {
				cout << "Level " << index + 1 << " from left to right: ";
				intuitionPrint(src, index);
				cout << endl;
			}
			else {
				cout << "Level " << index + 1 << " from right to left: ";
				reversePrint(src, index);
				cout << endl;
			}
		}
	}
	void printOPT(void) {
		queue<nod*> src;
		nod *tmp = nullptr;
		vector<nod*> print;
		vector<int> level_nodes(2, 0);
		int level = 0;

		src.push(root);
		level_nodes[level]++;
		while (!src.empty()) {
			tmp = src.front();
			src.pop();
			print.push_back(tmp);
			if (tmp->lchild) {
				src.push(tmp->lchild);
				level_nodes[level + 1]++;
			}
			if (tmp->rchild) {
				src.push(tmp->rchild);
				level_nodes[level + 1]++;
			}

			if (--level_nodes[level] == 0) {
				level++;
				level_nodes.push_back(0);

				if (level % 2) {
					cout << "Level " << level << " from left to right: ";
					for (size_t index = 0; index != print.size(); ++index) {
						cout << print[index]->value << ' ';
					}
					cout << endl;
				}
				else {
					cout << "Level " << level << " from right to left: ";
					for (int index = print.size() - 1; index >= 0; --index) {
						cout << print[index]->value << ' ';
					}
					cout <<  endl;
				}
				print.clear();
			}
		}
	}

/****************************test code****************************/
/*
	binaryTree x;
	vector<int> src{ -3, 3, -9, 1, 0, 2, 1, '#', '#', 1, 6, '#', '#', '#', '#', '#', '#', '#', '#' };
	node *root = x.generateFixedBinaryTree(src);
	x.binaryTreePrint();
	zigzagPrint tree(root);
	tree.print();
****debug:    -------------------
			  |                 |
			--3----           --'--
			|     |           |   |
			1   --0--         2   1
				|   |
				1   6
			Level 1 from left to right: -3
			Level 2 from right to left: -9 3
			Level 3 from left to right: 1 0 2 1
			Level 4 from right to left: 6 1
			nice!
*/
/*****************************end test****************************/
};
//#32
class modifySearchTree {
private:
	nod *src_root;
	nod *res_root;
	vector<int> src;
	binaryTree res_tree;
public:
	modifySearchTree(nod *_src_root = nullptr) : src_root(_src_root), res_root(nullptr) {}
	vector<nod*> getChanges(void) {
		nod *cur = src_root, *most_right = nullptr, *tmp = new nod(INT_MIN);
		vector< vector<nod*> > res(1);
		size_t index_res = 0;

		while (cur != nullptr) {
			if (cur->lchild) {
				most_right = cur->lchild;
				while (most_right->lchild != nullptr && most_right->lchild != cur)
					most_right = most_right->rchild;
				if (most_right->rchild == nullptr) {
					most_right->rchild = cur;
					cur = cur->lchild;
					continue;
				}
				else {
					most_right->rchild = nullptr;
					src.push_back(cur->value);
					if (tmp->value > cur->value) {
						res[index_res].push_back(tmp);
						res[index_res++].push_back(cur);
						res.push_back({});
					}
					tmp = cur;
				}
			}
			else {
				src.push_back(cur->value);
				if (tmp->value > cur->value) {
					res[index_res].push_back(tmp);
					res[index_res++].push_back(cur);
					res.push_back({});
				}
				tmp = cur;
			}
			cur = cur->rchild;
		}
		tmp = nullptr;
		delete(tmp);
		return{ res[0][0] , res[res.size() - 2][1] };
	}
/****************************test code****************************/
/*
    binaryTree x;
	vector<int> src1{ 6, 4, 10, 2, 3, 9, 13,  '#', '#', '#', '#', 7, 8, '#', '#', '#', '#', '#', '#', '#', '#' };
	vector<int> src2{ 18, 6, 15, 3, 7, 12, 10,  '#', '#', '#', '#', 11, 14, '#', '#', '#', '#', '#', '#', '#', '#', '#', '#' };
	node *root = x.generateFixedBinaryTree(src2);
	x.binaryTreePrint();
	modifySearchTree tree(root);
	vector<node*> res = tree.getChanges();
	cout << res[0]->value << ' ' << res[1]->value << endl;
	//tree.modify();
****debug:    ----18-----------
			  |               |
			--6--        -----15--
			|   |        |       |
			3   7     ---12--    10
					  |     |
					  11    14
			18 10
			nice!
*/
/*****************************end test****************************/
};
//#33  
class compareTopoTree {
private:
	nod *src_root1;
	nod *src_root2;
private:
	bool check(nod *_node, nod *_root) {
		if (_root == nullptr)
			return true;
		if ((_node == nullptr && _root != nullptr) || _node->value != _root->value)
			return false;
		return check(_node->lchild, _root->lchild) && check(_node->rchild, _root->rchild);
	}
	nod *find(nod *_head) {
		queue<nod*> que;
		nod *tmp = nullptr;

		que.push(_head);
		while (!que.empty()) {
			tmp = que.front();
			if (tmp->value == src_root2->value)
				return tmp;
			que.pop();
			if (tmp->lchild)
				que.push(tmp->lchild);
			if (tmp->rchild)
				que.push(tmp->rchild);
		}
	}
public:
	compareTopoTree(nod *_src_root1 = nullptr, nod *_src_root2 = nullptr) :
		src_root1(_src_root1), src_root2(_src_root2) {}
	bool getResult(void) {
		queue<nod*> que;
		nod *find_res = nullptr, *tmp = nullptr;

		que.push(src_root1);
		while (!que.empty()) {
			tmp = que.front();
			find_res = find(tmp);
			if (check(find_res, src_root2) == true)
				return true;
			else {
				que.pop();
				if (tmp->lchild)
					que.push(tmp->lchild);
				if (tmp->rchild)
					que.push(tmp->rchild);
			}
		}
		if (que.empty())
			return false;
	}
/****************************test code****************************/
/*
	binaryTree x;
	vector<int> src1{ 1, 2, 3, 4, 5, 6, 7, '#', 9, 10, '#', 2, '#', '#', '#', '#', '#', '#', '#', 4, 5, 9, '#', '#', '#', '#', '#' };
	vector<int> src2{ 2, 4, 5, 8, '#', '#', '#', '#', '#'};
	node *root1 = x.generateFixedBinaryTree(src1);
	x.binaryTreePrint();
	cout << endl << endl;
	node *root2 = x.generateFixedBinaryTree(src2);
	x.binaryTreePrint();
	compareTree c(root1, root2);
	cout << (c.getResult() ? "yeah!" : "nope") << endl;
****debug:      -------1-------------------
				|                         |
			----2-----                  --3--
			|        |                  |   |
			4--   ---5              ----6   7
			  |   |                 |
			  9   10              --2--
								  |   |
								--4   5
								|
								9


			  --2--
			  |   |
			--4   5
			|
			8
			nope
			nice!
*/
/*****************************end test****************************/
};
//#34
class compareTopoSubTree {
private:
	nod *src_root;
	nod *cmp_root;
	string src_serial;
	string cmp_serial;
private:
	string preSerial(nod *_head) {
		if (_head == nullptr)
			return "#!";
		string res;
		res.push_back(_head->value + '0');
		res += '!';
		res += preSerial(_head->lchild);
		res += preSerial(_head->rchild);
		return res;
	}
	vector<int> getNextArray(const string _tmp) {

	}
	int getIndefOf(const string &_src, const string &_cmp) {
		size_t res = _src.find(_cmp);
		return res;
	}
public:
	compareTopoSubTree(nod *_src_root = nullptr, nod *_cmp_root = nullptr) :
		src_root(_src_root), cmp_root(_cmp_root) {}
	bool getResult(void) {
		src_serial = preSerial(src_root);
		cmp_serial = preSerial(cmp_root);
		return (src_serial.find(cmp_serial) == string::npos ? false : true);
	}





/****************************test code****************************/
/*
    binaryTree x;
	vector<int> src1{ 1, 2, 3, 4, 5, 6, 7, '#', 9, 10, '#', 2, '#', '#', '#', '#', '#', '#', '#', 4, 5, 8, '#', '#', '#', '#', '#' };
	vector<int> src2{ 2, 4, 5, 8, '#', '#', '#', '#', '#'};
	vector<int> src3{ 2, 4, 5, 8, 9, '#', '#', '#', '#', '#', '#' };
	node *root1 = x.generateFixedBinaryTree(src1);
	x.binaryTreePrint();
	cout << endl << endl;
	node *root2 = x.generateFixedBinaryTree(src2);
	x.binaryTreePrint();
	node *root3 = x.generateFixedBinaryTree(src3);
	x.binaryTreePrint();
	compareTopoSubTree c1(root1, root2);
	cout << (c1.getResult() == true ? "yeah!" : "nope") << endl;
	compareTopoSubTree c2(root1, root3);
	cout << (c2.getResult() == true ? "yeah!" : "nope") << endl;
****debug:      -------1-------------------
				|                         |
			----2-----                  --3--
			|        |                  |   |
			4--   ---5              ----6   7
			  |   |                 |
			  9   10              --2--
								  |   |
								--4   5
								|
								8


			  --2--
			  |   |
			--4   5
			|
			8
			  ----2----
			  |       |
			--4--     5
			|   |
			8   9
			yeah!
			nope
			nice!
*/
/*****************************end test****************************/
};
//#35
class isBlanceTree {
private:
	nod *root;
private:
	int getHeight(nod *_head) {
		if (_head == nullptr)
			return 0;
		getHeight(_head->lchild);
		getHeight(_head->rchild);
	}
public:
	isBlanceTree(nod *_root = nullptr) : root(_root) {}
	bool getResult(void) {
		int l_height, r_height = 0;
		if(root->lchild)
            l_height = getHeight(root->lchild);
		if (root->rchild)
			r_height = getHeight(root->rchild);
		return (abs(l_height - r_height) > 1 ? false : true);
	}




/****************************test code****************************/
/*

****debug:
*/
/*****************************end test****************************/
};
//#36
class unicode {
private:
	string src;
	unordered_map<char, int> res;
public:
	unicode(const string &_src) : src(_src) {}
	string getResult(void) {
		char tmp = src[0];
		int cnt = 0;
		string res;

		for (size_t index = 0; index < src.size(); ) {
			while (tmp == src[index++])
				++cnt;
			res.push_back(cnt + '0');
			res.push_back(tmp);
			cnt = 0;
			tmp = src[--index];
		}

		return res;
	}
};
//#37
class printCommonList {
private:
	list<int> src_head;
	list<int> cmp_head;
public:
	printCommonList(list<int> _src_head, list<int> _cmp_head) : src_head(_src_head), cmp_head(_cmp_head) {}
	void getResult(void) {
		list<int> res_head;

		cout << "src_head: ";
		listPrint(src_head);
		cout << "cmp_head: ";
		listPrint(cmp_head);
		while (!src_head.empty() && !cmp_head.empty()) {
			if (src_head.front() < cmp_head.front())
				src_head.pop_front();
			else if (src_head.front() > cmp_head.front())
				cmp_head.pop_front();
			else {
				res_head.push_back(cmp_head.front());
				src_head.pop_front();
				cmp_head.pop_front();
			}
		}
		cout << "res_head: ";
		listPrint(res_head);
	}
/****************************test code****************************/
/*
    list<int> src_head = { 1, 2, 3, 4, 5, 6, 7, 8 };
	list<int> cmp_head = {4, 5, 6, 7, 8, 9, 10};
	printCommonList p(src_head, cmp_head);
	p.getResult();
****debug:  src_head: 1 2 3 4 5 6 7 8
			cmp_head: 4 5 6 7 8 9 10
			res_head: 4 5 6 7 8
			nice!
*/
/*****************************end test****************************/
};
//#38
class bubbleSort {
private:
	vector<int> src;
public:
	bubbleSort(const vector<int> &_src) : src(_src) {}
	~bubbleSort(void) {}
	void sort(void) {
		for (int i = src.size() - 1; i > 0; --i) {
			for (int j = 0; j < i; ++j) {
				if (src[j] > src[j + 1])
					swap(src[j], src[j + 1]);
			}
		}
		vectorPrint<int>(src);
	}

/****************************test code****************************/
/*
	vector<int> src = {6, 3, 5, 1, 8, 2, 12, 4, 7, 14};
	vectorPrint<int>(src);

	bubbleSort x(src);
	x.sort();
****debug: 6 3 5 1 8 2 12 4 7 14
		   1 2 3 4 5 6 7 8 12 14
*/
/*****************************end test****************************/
};
//#39
class selectSort {
private:
	vector<int> src;
public:
	selectSort(const vector<int> _src) : src(_src) {}
	~selectSort(void) {}
	void sort(void) {
		int minIndex = 0;
		for (int i = 0; i < src.size() - 1; ++i) {
			minIndex = i;
			for (int j = i+1; j < src.size() - 1; ++j)
				minIndex = src[minIndex] < src[j] ? minIndex : j;
			swap(src[i], src[minIndex]);
		}
		vectorPrint<int>(src);
	}
/****************************test code****************************/
/*
	vector<int> src = { 6, 3, 5, 1, 8, 2, 12, 4, 7, 14 };
	vectorPrint<int>(src);

	selectSort x(src);
	x.sort();
****debug: 6 3 5 1 8 2 12 4 7 14
		   1 2 3 4 5 6 7 8 12 14
		   nice!
*/
/*****************************end test****************************/
};
//#40
class insertSort {
private:
	vector<int> src;
public:
	insertSort(const vector<int> _src) : src(_src) {}
	~insertSort(void) {}
	void sort(void) {
		for (int i = 1; i < src.size() - 1; ++i) {
			for (int j = i - 1; (j >= 0) && (src[j + 1] < src[j]); --j)
				swap(src[j], src[j + 1]);
		}
		vectorPrint<int>(src);
	}



/****************************test code****************************/
/*
	vector<int> src = { 6, 3, 5, 1, 8, 2, 12, 4, 7, 14 };
	vectorPrint<int>(src);

	insertSort x(src);
	x.sort();
****debug: 6 3 5 1 8 2 12 4 7 14
		   1 2 3 4 5 6 7 8 12 14
		   nice!
*/
/*****************************end test****************************/
};
//#41
class reverseList {
private:
	struct nod {
		int data;
		nod *next;
		nod(int _data = 0, nod *_next = nullptr) : data(_data), next(_next) {}
	};
	nod *head;
	vector<int> data;
public:
	reverseList(vector<int> _data) : data(_data), head(nullptr) {}
	~reverseList(void) {}
	void genersteList(void) {
		nod *cur = nullptr, *tmp = nullptr;
		head = new nod(data[0]);
		tmp = head;
		for (size_t i = 1; i < data.size(); ++i) {
			cur = new nod(data[i]);
			tmp->next = cur;
			tmp = cur;
		}
		cur->next = nullptr;
	}
	void reverse1(void) {
		stack<nod*> nodes;
		nod *tmp = head, *cur = nullptr;
		while (tmp != nullptr) {
			nodes.push(tmp);
			tmp = tmp->next;
		}
		head = nodes.top();
		nodes.pop();
		head->next = nodes.top();
		while (!nodes.empty()) {
			tmp = nodes.top();
			nodes.pop();
			tmp->next = (nodes.empty() ? nullptr : nodes.top());
		}
	}
	void reverse2(void) {
		nod *front = nullptr, *cur = nullptr, *next = nullptr;

		cur = head;
		while (cur != nullptr) {
			next = cur->next;
			cur->next = front;
			front = cur;
			cur = next;
		}
		head = front;
	}


/****************************test code****************************/
/*
	vector<int> src = {1, 2, 3, 4, 5, 6, 7};
	reverseList x(src);

	x.genersteList();
	x.reverse1();
	x.reverse2();
****debug:
*/
/*****************************end test****************************/
};
//#42



//##
class myTest {
private:
	nod *root;
	vector<nod *>src;
public:
	myTest(nod *_root = nullptr) : root(_root) {}
	void morriss(void) {
		nod *cur = root, *most_right = nullptr;
		while (cur != nullptr) {
			src.push_back(cur);
			if (cur->lchild != nullptr) {
				most_right = cur->lchild;
				while (most_right->rchild != nullptr && most_right->rchild != cur)
					most_right = most_right->rchild;
				if (most_right->rchild == nullptr) {
					most_right->rchild = cur;
					cur = cur->lchild;
					continue;
				}
				else
					most_right->rchild = nullptr;
			}
			cur = cur->rchild;
		}
	}
	void preTraverse(void) {
		nod *cur = root, *most_right = nullptr;
		while (cur != nullptr) {
			if (cur->lchild != nullptr) {
				most_right = cur->lchild;
				while (most_right->rchild != nullptr && most_right->rchild != cur)
					most_right = most_right->rchild;
				if (most_right->rchild == nullptr) {
					most_right->rchild = cur;
					cout << cur->value << ' ';
					cur = cur->lchild;
					continue;
				}
				else
					most_right->rchild = nullptr;
			}
			else {
				cout << cur->value << ' ';
			}
			cur = cur->rchild;
		}
	}
	void midTraverse(void) {
		nod *cur = root, *most_right = nullptr;
		while (cur != nullptr) {
			if (cur->lchild != nullptr) {
				most_right = cur->lchild;
				while (most_right->rchild != nullptr && most_right->rchild != cur)
					most_right = most_right->rchild;
				if (most_right->rchild == nullptr) {
					most_right->rchild = cur;
					cur = cur->lchild;
					continue;
				}
				else {
					most_right->rchild = nullptr;
					cout << cur->value << ' ';
				}
			}
			else
				cout << cur->value << ' ';
			cur = cur->rchild;
		}
	}
	void postTraverse(void) {
		nod *cur = root, *most_right = nullptr;
		while (cur != nullptr) {
			if (cur->lchild != nullptr) {
				most_right = cur->lchild;
				while (most_right->rchild != nullptr && most_right->rchild != cur)
					most_right = most_right->rchild;
				if (most_right->rchild == nullptr) {
					most_right->rchild = cur;
					cur = cur->lchild;
					continue;
				}
				else
					most_right->rchild = nullptr;
			}
			cur = cur->rchild;
		}
	}
/****************************test code****************************/
/*
    binaryTree tree;
	node *root = tree.generateRandomBinaryTree(30, 30);
	tree.binaryTreePrint();
	myTest x(root);
	x.morriss();
	x.midTraverse();
	cout << endl;
	x.preTraverse();
	cout << endl;
	x.postTraverse();
****debug:        --------17-----------------
				  |                         |
			  ----4----                   --10-
			  |       |                   |   |
			--29-   --4--                 18  18
			|   |   |   |
			22  14  5   5
			22 29 14 4 5 4 5 17 18 10 18
			17 4 29 22 14 4 5 5 10 18 18
			nice!
*/
/*****************************end test****************************/
};



int main(void)
{
	int test_time = 1000;
	int max_size = 30;
	int max_value = 100;
	bool succeed = true;

	//binaryTree x;
	//vector<int> src1{ 1, 2, 3, 4, 5, 6, 7, '#', 9, 10, '#', 2, '#', '#', '#', '#', '#', '#', '#', 4, 5, 8, '#', '#', '#', '#', '#' };
	//vector<int> src2{ 2, 4, 5, 8, '#', '#', '#', '#', '#'};
	//vector<int> src3{ 2, 4, 5, 8, 9, '#', '#', '#', '#', '#', '#' };
	//node *root1 = x.generateFixedBinaryTree(src1);
	//x.binaryTreePrint();
	//cout << endl << endl;
	//node *root2 = x.generateFixedBinaryTree(src2);
	//x.binaryTreePrint();
	//node *root3 = x.generateFixedBinaryTree(src3);
	//x.binaryTreePrint();

	




	cout << (succeed == true ? "nice!" : "fucking fucked!!!") << endl;
	system("pause");
	return 0;
}