//---------------------------------------------------------------------------

#ifndef lista_H
#define lista_H
//---------------------------------------------------------------------------

struct node_double_linked{
  void* inf;
  node_double_linked* next, *prev;
};
//---------------------------------------------------------------------------
class TLista{
	private:
	public:
		node_double_linked* head, *tail;
		int count;
		TLista();
		TLista(TLista&);
		~TLista();
		void Add(void*);
		void Delete(int);
	node_double_linked* DeleteCurrent(node_double_linked*);
	void* GetInfo(int);
	node_double_linked* GetNode(int Index);
	void Append(TLista & source);
	void Insert(long, void*);
	void Clear(void);
	void DeleteHead(void);


	void* GetCurrentInfo(node_double_linked* p);
	void* GetHeadInfo(void);
	void* GetTailInfo(void);

	void* GetNextInfo(node_double_linked*);
	void* GetPrevInfo(node_double_linked*);
	void* GetNextCircularInfo(node_double_linked*);
	void* GetPrevCircularInfo(node_double_linked*);
	node_double_linked* delete_current_circular(node_double_linked* p);

	void AppendWithoutCopy(TLista &source);
	void make_circular(void);

};
//---------------------------------------------------------------------------
#endif
