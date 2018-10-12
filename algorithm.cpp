#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <stack>
#include <queue>
#include <list>
#include <unordered_map>
//#include <hash_map>

#include <string>
#include <windows.h>
#include <iomanip>

using namespace std;
using namespace stdext;

#define QUICKSORT     0
#define HEAPSORT      0
#define COMPARATOR    0
#define BUCKETSORT    0
#define PALINDROME    0
#define BINARYTREE    0
#define THREADED      0
#define HASHMAP       0
#define UNIONFINDSET  0
#define TRIETREE      0
#define MINIPATH      1

#define PRINTSPACE(n); for(int i=0; i<n; ++i) {cout << ' ';}
#define PRINTVALUE(space, value); { PRINTSPACE(space); cout << value;}
#define SERIAL(value); {serial += to_string(value);; serial.push_back(' ');}
//template<class T>
//static void comparator(vector<T> &arr);

template<class T>
static vector<T> generateRandomArray(int maxsize, int maxvalue)
{
	int size = (int)((maxsize + 1) * (rand() % 10 / 10.0));
	vector<T> arr;
	for (int i = 0; i < size; ++i)
		arr.push_back((T)((maxvalue + 1) * (rand() % 10 / 10.0)));// -(T)(maxvalue * (rand() % 10 / 10.0)));

	return arr;
}

template<class T>
static bool isEqual(const vector<T> &arr1, const vector<T> &arr2)
{
	if (arr1.size() == 0 || arr2.size() == 0)
		return false;
	if (arr1.size() != arr2.size())
		return false;
	for (size_t i = 0; i < arr1.size(); ++i)
	{
		if (arr1[i] != arr2[i])
			return false;
	}

	return true;
}

template<class T>
static vector<int> copyArray(const vector<T> &arr)
{
	vector<T> tmp;

	if (arr.size() == 0)
		return tmp;
	for (size_t i = 0; i < arr.size(); ++i)
		tmp.push_back(arr[i]);

	return tmp;
}

template<class T>
static void printArray(const vector<T> &arr)
{
	if (arr.size() == 0)
		return;
	for (size_t i = 0; i < arr.size(); ++i)
		cout << arr[i] << ' ';
	cout << endl;
}
/****************************test code******************************/
#if QUICKSORT
template<class T>
static void quickSort(vector<T> &arr);
template<class T>
static void quickSort(vector<T> &arr, int l, int r);
template<class T>
static vector<T> parition(vector<T> &arr, int l, int r);

template<class T>
static void comparator(vector<T> &arr) {
	sort(arr.begin(), arr.end());
}

template<class T>
static void quickSort(vector<T> &arr)
{
	if (arr.size() < 2)
		return;
	quickSort(arr, 0, arr.size() - 1);
}

template<class T>
static void quickSort(vector<T> &arr, int l, int r)
{
	if (l < r)
	{
		vector<T> tmp = parition(arr, l, r);
		quickSort(arr, l, tmp[0] - 1);
		quickSort(arr, tmp[1] + 1, r);
	}
}

template<class T>
static vector<T> parition(vector<T> &arr, int l, int r)
{
	int less = l - 1;
	int more = r + 1;
	int num = arr[r];
	int cur = l;
	vector<T> tmp;

	while (cur < more)
	{
		if (arr[cur] < num)
			swap(arr[++less], arr[cur++]);
		else if (arr[cur] > num)
			swap(arr[--more], arr[cur]);
		else
			++cur;
	}
	tmp.push_back(less + 1);
	tmp.push_back(more - 1);

	return tmp;
}
#elif HEAPSORT
template<class T>
static void heapInsert(vector<T> &arr, int i);
template<class T>
static void heapify(vector<T> &arr, int index, int size);

template<class T>
static void comparator(vector<T> &arr) {
	sort(arr.begin(), arr.end());
}

template<class T>
static void heapSort(vector<T> &arr) {
	if (arr.size() < 2)
		return;
	for (size_t i = 0; i < arr.size(); ++i)
		heapInsert(arr, i);
	size_t size = arr.size();
	swap(arr[0], arr[--size]);
	while (size > 0) {
		heapify(arr, 0, size);
		swap(arr[0], arr[--size]);
	}
}

template<class T>
static void heapInsert(vector<T> &arr, int index) {
	while (arr[index] > arr[(index - 1) / 2]) {
		swap(arr[index], arr[(index - 1) >> 1]);
		index = ((index - 1) >> 1);
	}
}

template<class T>
static void heapify(vector<T> &arr, int index, int size) {
	int left = (index >> 1) + 1;
	while (left < size) {
		int largest = (left + 1 < size) && arr[left + 1] > arr[left] ? left + 1 : left;
		largest = arr[largest] > arr[index] ? largest : index;
		if (largest == index)
			break;
		swap(arr[largest], arr[index]);
		index = largest;
		left = (index << 1) + 1;
	}
}

template<class T>
static void comparator(vector<T> &arr) {
	sort(arr.begin(), arr.end());
}
#elif COMPARATOR
template<class T>
static void comparator(vector<T> &arr) {
	sort(arr.begin(), arr.end());
}
class Student {
public:
	string name;
	int id;
	int age;
public:
	Student(string _name, int _id, int _age) : name(_name), id(_id), age(_age) {}
};
static bool IDComparator(Student &s1, Student &s2) {
	return (s1.id < s2.id);
}
static bool AGEComparator(Student &s1, Student &s2) {
	return (s1.age < s2.age);
}
#elif BUCKETSORT
template<class T>
static void comparator(vector<T> &arr) {
	sort(arr.begin(), arr.end());
}
template<class T>
static void bucketSort(vector<T> &arr) {
	if (arr.size() < 2)
		return;
	int size = MININT;
	for (size_t i = 0; i < arr.size(); ++i)
		size = max(size, arr[i]);
	vector<T> bucket(size + 1);
	for (size_t i = 0; i < arr.size(); ++i)
		bucket[arr[i]]++;
	for (size_t i = 0, j = 0; j < size + 1; ++j) {
		while (bucket[j]--)
			arr[i++] = j;
	}
}
#elif PALINDROME
class node {
public:
	int value;
	node *next;
	node(int data, node *next = nullptr) { value = data; }
};
static bool isPalindromeN(node &head) {
	stack<node> tmp;
	node *cur = &head;
	while (cur != nullptr) {
		tmp.push(*cur);
		cur = cur->next;
	}

	cur = &head;
	while (cur != nullptr) {
		node x = tmp.top();
		tmp.pop();
		if (cur->value != x.value)
			return FALSE;
		cur = cur->next;
	} 
	return TRUE;
}

static bool isPalindromeN2(node &head) {
	node *r2l = nullptr,  *l2r = nullptr;
	node *size = nullptr, *tmp = nullptr, *mid = nullptr;

	l2r  = &head;
	size = &head;
	mid  = &head;
	while ((size->next != nullptr) && (size->next->next != nullptr) && (mid->next != nullptr)) {
		size = size->next->next;
		mid = mid->next;
	}
	mid->next = nullptr;
	l2r = head.next;
	(*tmp).next = nullptr;
	while (l2r != mid) {
		tmp = l2r->next;
		(*tmp).next = l2r;
		l2r = l2r->next;
	}
	return TRUE;
}
#elif BINARYTREE
class BiTree {
private:
	struct node {
		int value;
		node *left, *right;
		node() : left(nullptr), right(nullptr) {}
		node(int _value, node *_l = nullptr, node *_r = nullptr) : value(_value), left(_l), right(_r) {  }
	};
	node *root;
	string serial;
	int depth;
	queue<node *> src;
public:
	BiTree(node *_root = nullptr, int _depth = 0) : root(_root) {}
	void getDepth(void);
	void creatBinary(void);
	void preOrderTraverse(void);
	void inOrderTraverse(void);
	void postOrderTraverse(void);
	void printBinary(void);
	void serialBinary(void);
};

void BiTree::getDepth(void) {
	if (root == nullptr) {
		return;
	}
	queue<node *> q;
	q.push(root);
	while (!q.empty()) {
		int len = q.size();
		depth++;
		while (len--) {
			node* tmp = q.front();
			q.pop();
			if (tmp->left) 
				q.push(tmp->left);
			if (tmp->right) 
				q.push(tmp->right);
		}
	}
}

void BiTree::creatBinary(void) {
	node *tmp = nullptr;
	int lchild = 0, rchild = 0, x = 0;
	queue<node *> que;
	int reserve = 0;

	cout << "please enter root value: ";
	cin >> x;
	if (x != -1) {
		root = new node(x);
		que.push(root);
		src.push(root);
		SERIAL(root->value);
		while (!que.empty()) {
			tmp = que.front();
			que.pop();
			cout << "please enter " << tmp->value << "'s lchild and rchild: ";
			cin >> lchild >> rchild;
			if (lchild != -1) {
				que.push(tmp->left = new node(lchild));
				src.push(new node(lchild));
				SERIAL(lchild);
			}
			else {
				src.push(new node('#'));
				serial += "# ";
			}
			if (rchild != -1) {
				que.push(tmp->right = new node(rchild));
				src.push(new node(rchild));
				serial += "# ";
			}
			else {
				src.push(new node('#'));
				serial += "#";
			}
		}
	}
}

void BiTree::preOrderTraverse(void) {
	stack<node *> tmp;
	node *x = root;

	cout << "preOrderTraverse:\n";
	tmp.push(x);
	while (!tmp.empty()) {
		x = tmp.top();
		tmp.pop();
		cout << x->value << ' ';
		if (x->right)
			tmp.push(x->right);
		if (x->left)
			tmp.push(x->left);
	}
	cout << endl;
}

void BiTree::inOrderTraverse(void) {
	stack<node *> tmp;
	node *x = nullptr;

	cout << "inOrderTraverse:\n";
	tmp.push(root);
	while (!tmp.empty()) {
		while ((x = tmp.top()) && x != nullptr)
			tmp.push(x->left);
		tmp.pop();
		if (!tmp.empty()) {
			x = tmp.top();
			tmp.pop();
			cout << x->value << ' ';
			tmp.push(x->right);
		}
	}
	cout << endl;
}
 
void BiTree::postOrderTraverse(void) {
	stack<node *> tmp;
	node *x = nullptr;

	cout << "postOrderTraverse:\n";
	tmp.push(root);
	x = tmp.top();
	while ((x->left) || x->right){
		if (x->right)
			tmp.push(x->right);
		if (x->left)
			tmp.push(x->left);
		x = tmp.top();
	}
	while (!tmp.empty()) {
		x = tmp.top();
		cout << x->value << ' ';
		tmp.pop();
	}
	cout << endl;
}

void BiTree::printBinary(void) {
	node *tmp = root;

	cout << "printBT:\n\n\n";
	/*get level*/
	this->getDepth();
	/*print nodes*/
	int position = ((1 << depth) >> 1);
	tmp = src.front();
	src.pop();
	PRINTSPACE(position);
	cout << tmp->value << endl;
	for (int j = 2; j < depth + 1; ++j) {
		position = 1 << (depth - j);
		for (int i = 0; i < (1 << (j - 1)); ++i) {
			tmp = src.front();
		    src.pop();
			PRINTSPACE(position);
			if (tmp->value != '#')
				cout << tmp->value;
			else
				cout << ' ';
			position = 1 << (depth + 1 - j);
			--position;
		}
		cout << endl;
	}
}

void BiTree::serialBinary(void) {
	cout << "serialBinary...\n";
	cout << serial << endl;
}

#elif THREADED
class threadedBT {
private:
	struct node {
		int value;
		int ltag;
		int rtag;
		node *left;
		node *right;
	public:
		node(int _val, int _ltag = 0, int _rtag = 0, node *_l = nullptr, node *_r = nullptr)
			: value(_val), left(_l), right(_r), ltag(_ltag), rtag(_rtag) {}
		node *getPredcessorNode(void) { return this->left; }
		node *getSuccessorNode(void)  { return this->right; }
	};
	node *root;
	node *head;
	int  levels;
public:
	threadedBT() : root(nullptr), head(nullptr), levels(0) {}
	void creatBT(void);
	void threadBT(void);
	void traverseBT(void);
};
void threadedBT::creatBT(void) {
	queue<node *> q;
	node *tmp = nullptr;
	int x = 0, ldata = 0, rdata = 0;

	head = new node(-1);
	cout << "createBT!\n";
	cout << "please enter root node: ";
	cin >> x;
	if (x != -1) {
		++levels;
		root = new node(x);
		q.push(root);
		while (!q.empty()) {
			tmp = q.front();
			q.pop();
			cout << "please enter " << tmp->value << "'s ldata and rdata: ";
			cin >> ldata >> rdata;
			if (ldata != -1 || rdata != -1)
				++levels;
			if (ldata != -1)
				q.push(tmp->left = new node(ldata));
			else
				tmp->ltag = 1;
			if (rdata != -1)
				q.push(tmp->right = new node(rdata));
			else
				tmp->rtag = 1;
		}
	}
}
void threadedBT::threadBT(void) {
	stack<node *> s;
	stack<node *> x;
	node *tmp  = nullptr;
	node *bkp  = head;

	cout << "threadedBT...\n";
	s.push(root); 
	while (!s.empty()) {
		while ((tmp = s.top()) && tmp)
			s.push(tmp->left);
		s.pop();
		if (!s.empty()) {
			tmp = s.top();
			s.pop();
			s.push(tmp->right);
			if (bkp->rtag == 1 )
				bkp->right = tmp;
			if (tmp->ltag == 1 )
			    tmp->left  = bkp;
			bkp = tmp;
		}
	}
	bkp->right = head;
}

void threadedBT::traverseBT(void) {
	node *tmp = root;
	
	cout << "traverseBT:";
	while (tmp != head) {
		while (tmp->ltag == 0)
			tmp = tmp->left;
		cout << tmp->value << ' ';
		while (tmp->rtag == 1 && tmp->right != head) {
			tmp = tmp->right;
			cout << tmp->value << ' ';
		}
		tmp = tmp->right;
	}
	cout << endl;
}
#elif HASHMAP
struct str_hash {
	size_t operator()(const string& str) const
	{
		unsigned long __h = 0;
		for (size_t i = 0; i < str.size(); i++)
			__h = 5 * __h + str[i];
		return size_t(__h);
	}
	bool operator()(const string& p1, const string& p2) const {
		return p1 == p2;
	}
};
#elif UNIONFINDSET
struct node {
	int value;
	node(int _value = 0) : value(_value) {}
	bool operator!=(const node &_node) {
		return (this->value != _node.value);
	}
	bool operator==(const node &_node) {
		return (this->value == _node.value);
	}
};
unordered_map<node, node> fatherMap;  //key: child    value: father
unordered_map<node, int> sizeMap;     //key: node     value: node's set node number
node findHead(node _node) {
	node father = fatherMap.find(_node)->first;
	if (father != _node) {
		father = fatherMap.find(father)->first;
	}
	fatherMap.insert(pair<node, node>(_node, father));    //optimization
	return father;
}
static void makeSet(list<node> &_nodes) {
	fatherMap.clear();
	sizeMap.clear();
	for (node x : _nodes) {
		fatherMap.insert(pair<node, node>(x, x));
		sizeMap.insert(pair<node, int>(x, 1));
	}
}
static bool isSameSet(const node &_a, const node &_b) {
	return (findHead(_a) == findHead(_b));
}
static void unionSet(const node &_a, const node &_b) {
	node head_a = findHead(_a);
	node head_b = findHead(_b);
	if (head_a == head_b)
		return;
	int size_a = sizeMap.find(head_a)->second;
	int size_b = sizeMap.find(head_b)->second;
	if (size_a <= size_b) {
		fatherMap.insert(pair<node, node>(head_a, head_b));
		sizeMap.insert(pair<node, int>(head_b, size_a + size_b));
	}	
	else {
		fatherMap.insert(pair<node, node>(head_b, head_a));
		sizeMap.insert(pair<node, int>(head_a, size_a + size_b));
	}
}

#elif TRIETREE
class trieTree {
private:
	struct node {
		int path;
		int end;
		node *next[26];
		node(int _path = 0, int _end = 0) : path(_path), end(_end) { memset(next, 0, sizeof(next)); }
	};
	node root;
private:
	void clearNode(node &_tmp) {
		for (int i = 0; i < 26; ++i) {
			if (_tmp.next[i] == nullptr)
				continue;
			clearNode(*(_tmp.next[i]));
			delete _tmp.next[i];
		}
	}
public:
	trieTree(void) {}
	~trieTree(void) { clearNode(root); }
	void insertStr(const string &_word) {
		size_t i = 0;
		int index = 0;
		node *tmp = &root;

		for (i = 0; i < _word.length(); ++i) {
			index = _word[i] - 'a';
			if (tmp->next[index] == nullptr) {
				tmp->next[index] = new node();
			}
			tmp = tmp->next[index];
			tmp->path++;
		}
		tmp->end++;
	}
	void deleteStr(const string &_word) {
		size_t i = 0;
		int index = 0;
		node *tmp = &root;

		if (searchStr(_word) == false)
			return;
		for (i = 0; i < _word.length(); ++i) {
			index = _word[i] - 'a';
			if (tmp->next[index]->path-- == 0)
				delete tmp->next[index];
				tmp = tmp->next[index];
		}
		tmp->end--;
	}
	bool searchStr(const string &_word) {
		size_t i = 0;
		int index = 0;
		node *tmp = &root;

		for (i = 0; i < _word.length(); ++i) {
			index = _word[i] - 'a';
			if (tmp->next[index]) {
				tmp = tmp->next[index];
				continue;
			}
			return false;
		}
		return tmp->end ? true : false;
 	}
	int preFixNumber(const string &_word) {
		size_t i = 0;
		int index = 0;
		node *tmp = &root;

		for (i = 0; i < _word.length(); ++i) {
			index = _word[i] - 'a';
			if (tmp->next[index]) {
				tmp = tmp->next[index];
			}
			else
			    return 0;
		}
		return tmp->path;
	}
/*************************test*********************/
	/*
	trieTree t;

	t.insertStr("abc");
	t.insertStr("cd");
	t.insertStr("cde");
	t.insertStr("cdef");
	cout << (t.searchStr("cde") ? "src has <cde>\n" : "src don't have cde\n");
	t.deleteStr("cde");
	cout << (t.searchStr("cde") ? "src has <cde>\n" : "src don't have cd\n");

	cout << "preFixNumber <cd>: " << t.preFixNumber("cd") << endl;
	debug: src has <cde>
	src don't have cd
	preFixNumber <cd>: 2
*/
};
#elif MINIPATH
class miniPath {
private:
	int(*src)[4];
	int col, row;
	int path;
private:
	int getminiPathRec(int(*_src)[4], int _row, int _col) {
		int tmp = _src[_row][_col];
		if ((_row == row - 1) && (_col == col - 1))
			return tmp;
		if (_row == row - 1)  //last row
			return tmp + getminiPathRec(_src, _row, _col + 1);
		if (_col == col - 1)  //last column
			return tmp + getminiPathRec(_src, _row + 1, _col);
		int right = getminiPathRec(_src, _row, _col + 1);
		int down  = getminiPathRec(_src, _row + 1, _col);
		return ((right < down ? right : down) + tmp);
	}
	int getminiPathDP(int(*src)[4]) {

	}
public:
	miniPath(int(*_src)[4] = nullptr, int _col = 0, int _row = 0) 
		: col(_col), row(_row), src(_src), path(0) {}
	int getminiPathRec(void) {
		return getminiPathRec(src, 0, 0);
	}
	int getminiPathDP(void) {
		return getminiPathDP(src);
	}
};
#endif
/*******************************End test******************************/
int main(void)
{
#if 0
	int testTime = 10000;
	int maxSize = 100;
	int maxValue = 30;
	bool succeed = true;
	int num = 6;
	vector<int> arr1, arr2;

	long timeStart = GetTickCount();
	for (int i = 0; i < testTime; ++i)
	{
		arr1 = generateRandomArray<int>(maxSize, maxValue);
		if (arr1.size() == 0)
			continue;
		arr2 = copyArray(arr1);

		//comparator(arr1);
		//bucketSort(arr2);

		if (!isEqual(arr1, arr2))
		{
			succeed = false;
			printArray(arr1);
			printArray(arr2);
			break;
		}
	}
	long timeStop = GetTickCount();
	cout << (succeed ? "Nice!\n" : "Fucking fucked!\n") << "Time: " << (timeStop - timeStart) << "ms\n";
#endif	
	int src[4][4] = { {1, 3, 5, 9},
	                  {8, 1, 3, 4}, 
					  {5, 0, 6, 1},
					  {8, 8, 4, 0}, };
	miniPath m(src, 4, 4);

	cout << "getminiPathRec: " << m.getminiPathRec() << endl;


	system("pause");
	return 0;
}