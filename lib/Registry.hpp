#ifndef UTILITIES_REGISTRY_HPP
#define UTILITIES_REGISTRY_HPP

#include <map>

namespace Utilities {

template <class TKey, class TValue>
struct Registry
{
	 static std::map<TKey,TValue>& registry()
	 {
		 static std::map<TKey,TValue> s;
		 return s;
	 }

	 static const TValue& query(const TKey &key, const TValue &defaultValue)
	 {
		 auto it = registry().find(key);
		 if (it == registry().end())
			 return defaultValue;
		 else
			 return it->second;
	 }

	 static const TValue& get(const TKey &key)
	 {
		 return registry().at(key);
	 }

	 static void registerItem(const TKey &key, const TValue &value)
	 {
		 registry().insert(std::make_pair(key,value));
	 }
};

} // namespace Utilities

#if 1
#define REGISTER_ITEM(type,id,item) \
	static int _reg_item_##item() { \
	Utilities::Registry<std::string,type>::registerItem(id,item); \
	return 0; } \
	int _var_reg_item_##item = _reg_item_##item(); 
	//__pragma(comment(linker, "/include:" "_var_reg_item_" #item))
#else
#define REGISTER_ITEM2(type,id,item) \
	static int _reg_item ## id () { \
		Utilities::Registry<std::string,type>::registerItem(#id,item); \
		return 0; \
	} \
	static int _var_reg_item ## id = _reg_item ## id (); 
#endif

#endif // UTILITIES_REGISTRY_HPP
